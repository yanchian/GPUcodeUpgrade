/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "ghost.hpp"
#include <sys/times.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

bool halt = false;

void signalHandler(int signal = 0) {
  halt = true;
}


// Main function on CPU /////////////////////////////////////////////////////////
int main(int argc, char** argv) {

  // capture SIGINT
  signal(SIGINT, signalHandler);

  if (argc < 2) {
    std::cerr << "Invoke with " << argv[0] << " <configuration file>" << std::endl;
    exit(1);
  }

  // select a device
  int num_devices;
  cudaGetDeviceCount(&num_devices);
  std::cout << "#Found " << num_devices << " GPGPUs" << std::endl;
  cudaDeviceProp properties;
  int best_device = 0;
  if (num_devices > 1) {
    // if there's more than one, pick the one with the highest compute capability
    int best_computemode = 0, computemode;
    for (int device = 0; device < num_devices; device++) {
      cudaGetDeviceProperties(&properties, device);
      std::cout << "  #" << device << " " << properties.name << ": " << properties.multiProcessorCount << " processors, compute capability " << properties.major << "." << properties.minor << std::endl;
      computemode = properties.major << 4 + properties.minor;
      if (best_computemode < computemode) {
        best_computemode = computemode;
        best_device = device;
      }
    }
  }
        best_device = atoi(argv[2]);
  cudaGetDeviceProperties(&properties, best_device);
  std::cout << "#  using #" << best_device << " (" << properties.name << ")" << std::endl;
  cudaSetDevice(best_device);

  // start a timer to get the total wall time at the end of the run
  struct tms startTimes, endTimes;
  timespec startClock, endClock;
  times(&startTimes);
  clock_gettime(CLOCK_REALTIME, &startClock);

  copyMechanismToDevice();

  Solver solver(argv[1]);

  // ============================================================
  // 初始条件诊断 (一次性)
  // ============================================================
  {
    Mesh<CPU>::type uCPU_init(*solver.u, Mesh<CPU>::type::Allocate);
    uCPU_init = *solver.u;

    int ig_x = (int)(ex_x * cell_x);
    int ig_y = (int)(ex_y * cell_y);
    int amb_x = (int)(0.1 * cell_x);  // 环境区采样点
    int amb_y = (int)(0.1 * cell_y);

    auto calcTemp = [](const real* q) -> real {
      real sumYinvM = 0.0, sumY_t = 0.0;
      for (int k = 0; k < NUM_SPECIES - 1; k++) {
        real Yk = q[SPECIES_INDEX(k)];
        if (Yk > 0.0) sumYinvM += Yk / MOL_WT_HOST[k];
        sumY_t += Yk;
      }
      real Ylast = 1.0 - sumY_t;
      if (Ylast > 0.0) sumYinvM += Ylast / MOL_WT_HOST[NUM_SPECIES - 1];
      real R_mix = R_UNIV * sumYinvM;
      return q[PRESSURE] / (q[DENSITY] * R_mix + 1e-30);
    };

    real q_ig[NUMBER_VARIABLES], q_amb[NUMBER_VARIABLES];
    conservativeToPrimitive(uCPU_init(ig_x, ig_y), q_ig);
    conservativeToPrimitive(uCPU_init(amb_x, amb_y), q_amb);

    std::cout << "# === 初始条件验证 ===" << std::endl;
    std::cout << "# 点火区 (" << ig_x << "," << ig_y << "):"
              << " P=" << q_ig[PRESSURE] << " Pa"
              << " rho=" << q_ig[DENSITY] << " kg/m3"
              << " T=" << calcTemp(q_ig) << " K"
              << " e_int=" << q_ig[PRESSURE]/((GAMMA_DEFAULT-1.0)*(q_ig[DENSITY]+1e-30))
              << " J/kg" << std::endl;
    std::cout << "#   物种: Y_H=" << q_ig[SPECIES_INDEX(SP_H)]
              << " Y_O=" << q_ig[SPECIES_INDEX(SP_O)]
              << " Y_OH=" << q_ig[SPECIES_INDEX(SP_OH)]
              << " Y_CH4=" << q_ig[SPECIES_INDEX(SP_CH4)]
              << " Y_O2=" << q_ig[SPECIES_INDEX(SP_O2)]
              << " sumY=" << (q_ig[SPECIES_INDEX(SP_H)]+q_ig[SPECIES_INDEX(SP_O)]+q_ig[SPECIES_INDEX(SP_OH)]+q_ig[SPECIES_INDEX(SP_CH4)]+q_ig[SPECIES_INDEX(SP_O2)])
              << std::endl;
    std::cout << "#   守恒量: DENSITY=" << uCPU_init(ig_x, ig_y)[DENSITY]
              << " ENERGY=" << uCPU_init(ig_x, ig_y)[ENERGY]
              << " E/rho=" << uCPU_init(ig_x, ig_y)[ENERGY]/uCPU_init(ig_x, ig_y)[DENSITY]
              << std::endl;
    std::cout << "# 环境区 (" << amb_x << "," << amb_y << "):"
              << " P=" << q_amb[PRESSURE] << " Pa"
              << " rho=" << q_amb[DENSITY] << " kg/m3"
              << " T=" << calcTemp(q_amb) << " K"
              << " Y_CH4=" << q_amb[SPECIES_INDEX(SP_CH4)]
              << " Y_O2=" << q_amb[SPECIES_INDEX(SP_O2)] << std::endl;
    std::cout << "# 期望: 点火区 P=4e4 rho=0.05 T≈2712K, 环境区 P=2e4 rho≈0.23 T≈293K" << std::endl;
    std::cout << "# ======================" << std::endl;
  }

  Solver::Status status = Solver::OUTPUT;

  Mesh<CPU>::type uCPU(*solver.u, Mesh<CPU>::type::Allocate);
  
#ifdef GLOUTPUT
  OpenGLOutputter outputter(argc, argv, *solver.u);
  boost::thread outputterThread(boost::ref(outputter));
#endif

  // open the data file and output a header line
  std::stringstream filename;
  filename << solver.outputDirectory << "data";

  std::ofstream dataFile;
  dataFile.open(filename.str().c_str());

#ifdef GHOST
  dataFile << "#" << solver.u->time() << "\t\t";
  for (int i = 0; i < solver.geometries.size(); i++) {
    if (solver.geometries[i].rotating) {
      dataFile << "torque on level set " << i << "\t\t";
    }
  }
#endif
  for (int i = 0; i < solver.outputRadii.size(); i++) {
    dataFile << "P(r=" << solver.outputRadii[i] << ")\t\t";
    dataFile << "flux(r=" << solver.outputRadii[i] << ")\t\t";
  }
  dataFile << std::endl;
  
   do {
	   
/*#ifdef GLOUTPUT
    if (solver.getStepNumber() % 1 == 0) {
      outputter.dispatchDraw(*solver.u);
      //outputter.paused = true;
      while (outputter.gridToRender != NULL); // SPIN
    }
#endif
*/

////////////////////////////////////////////////////////////////


    if (status == Solver::OUTPUT) {
    
      if (true) {  
	  uCPU = *solver.u;
        std::stringstream filename0;
        filename0 << solver.outputDirectory << "sensor" << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".txt";

        std::ofstream outFile0;
        outFile0.open(filename0.str().c_str());

        outFile0.precision(8);
		//outFile1 << "i" << " " << "j" << " " << "p" << ' ' << "rho" << ' ' << "t" << " " <<  "u" << " " <<  "v" << " " <<  "LA0" << " " <<  "LA1" << std::endl;
		//for (int j = uCPU.activeNy()/5; j < uCPU.activeNy(); j += uCPU.activeNy()/5) {
		int j = ex_y * cell_y;	
    int i = ex_x * cell_x + r_bump * 8 * cell_x / 9; 
			  
            real q[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), q);
            // 计算温度: T = P/(ρ*R_mix), 理想气体定律
            real sumYinvM = 0.0, sumY_t = 0.0;
            for (int k = 0; k < NUM_SPECIES - 1; k++) {
              real Yk = q[SPECIES_INDEX(k)];
              if (Yk > 0.0) sumYinvM += Yk / MOL_WT_HOST[k];
              sumY_t += Yk;
            }
            real Ylast = 1.0 - sumY_t;
            if (Ylast > 0.0) sumYinvM += Ylast / MOL_WT_HOST[NUM_SPECIES - 1];
            real R_mix = R_UNIV * sumYinvM;
            real T_out = q[PRESSURE] / (q[DENSITY] * R_mix + 1e-30);
            outFile0 << std::fixed << i << " " << j << " " << q[PRESSURE] << ' ' << q[DENSITY] << ' ' << T_out << " " <<  q[XVELOCITY] << " " <<  q[YVELOCITY] << " " <<  q[SPECIES_INDEX(SP_CH4)] << " " <<  q[SPECIES_INDEX(SP_O2)] << std::endl;
          
        outFile0.close(); 
    

        std::stringstream filename1;
        filename1 << solver.outputDirectory << "RESULT" << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".txt";

        std::ofstream outFile1;
        outFile1.open(filename1.str().c_str());

        outFile1.precision(8);
		//outFile1 << "i" << " " << "j" << " " << "p" << ' ' << "rho" << ' ' << "t" << " " <<  "u" << " " <<  "v" << " " <<  "YCH4" << " " <<  "YO2" << std::endl;
		//for (int j = uCPU.activeNy()/5; j < uCPU.activeNy(); j += uCPU.activeNy()/5) {
		j = ex_y * cell_y;	
          for (int i = 0; i < uCPU.activeNx(); ++i) {
			  
            real q[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), q);
            // 计算温度: T = P/(ρ*R_mix), 理想气体定律
            real sumYinvM = 0.0, sumY_t = 0.0;
            for (int k = 0; k < NUM_SPECIES - 1; k++) {
              real Yk = q[SPECIES_INDEX(k)];
              if (Yk > 0.0) sumYinvM += Yk / MOL_WT_HOST[k];
              sumY_t += Yk;
            }
            real Ylast = 1.0 - sumY_t;
            if (Ylast > 0.0) sumYinvM += Ylast / MOL_WT_HOST[NUM_SPECIES - 1];
            real R_mix = R_UNIV * sumYinvM;
            real T_out = q[PRESSURE] / (q[DENSITY] * R_mix + 1e-30);
            outFile1 << std::fixed << i << " " << j << " " << q[PRESSURE] << ' ' << q[DENSITY] << ' ' << T_out << " " <<  q[XVELOCITY] << " " <<  q[YVELOCITY] << " " <<  q[SPECIES_INDEX(SP_CH4)] << " " <<  q[SPECIES_INDEX(SP_O2)] << std::endl;
          }
        outFile1.close(); 

      // ============================================================
      // 诊断输出: 采样点火区和全域极值
      // ============================================================
      {
        auto calcTemp = [](const real* q) -> real {
          real sumYinvM = 0.0, sumY_t = 0.0;
          for (int k = 0; k < NUM_SPECIES - 1; k++) {
            real Yk = q[SPECIES_INDEX(k)];
            if (Yk > 0.0) sumYinvM += Yk / MOL_WT_HOST[k];
            sumY_t += Yk;
          }
          real Ylast = 1.0 - sumY_t;
          if (Ylast > 0.0) sumYinvM += Ylast / MOL_WT_HOST[NUM_SPECIES - 1];
          real R_mix = R_UNIV * sumYinvM;
          return q[PRESSURE] / (q[DENSITY] * R_mix + 1e-30);
        };

        // 采样点火区中心 (单点)
        int ig_x = (int)(ex_x * cell_x);
        int ig_y = (int)(ex_y * cell_y);

        real q1[NUMBER_VARIABLES];
        conservativeToPrimitive(uCPU(ig_x, ig_y), q1);
        real T1 = calcTemp(q1);

        // 全域扫描 min/max
        real rho_min = 1e30, rho_max = -1e30;
        real P_min = 1e30, P_max = -1e30;
        real T_min = 1e30, T_max = -1e30;
        real e_int_min = 1e30, e_int_max = -1e30;
        int n_neg_rho = 0, n_neg_P = 0;

        for (int j = 0; j < uCPU.activeNy(); j++) {
          for (int i = 0; i < uCPU.activeNx(); i++) {
            real q[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), q);
            real rho = q[DENSITY], P = q[PRESSURE];
            if (rho < 0) n_neg_rho++;
            if (P < 0) n_neg_P++;
            if (rho < rho_min) rho_min = rho;
            if (rho > rho_max) rho_max = rho;
            if (P < P_min) P_min = P;
            if (P > P_max) P_max = P;
            real T = calcTemp(q);
            if (T < T_min) T_min = T;
            if (T > T_max) T_max = T;
            real e_kin = 0.5 * (q[XVELOCITY]*q[XVELOCITY] + q[YVELOCITY]*q[YVELOCITY]);
            real e_int = P / ((GAMMA_DEFAULT - 1.0) * (rho + 1e-30));
            if (e_int < e_int_min) e_int_min = e_int;
            if (e_int > e_int_max) e_int_max = e_int;
          }
        }

        std::cout << "# DIAG step=" << solver.getStepNumber()
                  << std::setprecision(6) << std::scientific
                  << " t=" << solver.u->time() << " dt=" << solver.getLastDt() << std::endl;
        std::cout << "#   点火区 ("
                  << ig_x << "," << ig_y << "): T=" << T1
                  << " P=" << q1[PRESSURE] << " rho=" << q1[DENSITY]
                  << " e_int=" << q1[PRESSURE]/((GAMMA_DEFAULT-1.0)*(q1[DENSITY]+1e-30))
                  << " Y_H=" << q1[SPECIES_INDEX(SP_H)]
                  << " Y_O=" << q1[SPECIES_INDEX(SP_O)]
                  << " Y_OH=" << q1[SPECIES_INDEX(SP_OH)]
                  << " Y_CH4=" << q1[SPECIES_INDEX(SP_CH4)]
                  << " Y_O2=" << q1[SPECIES_INDEX(SP_O2)] << std::endl;
        std::cout << "#   全域: rho=[" << std::scientific << rho_min << ", " << rho_max
                  << "] P=[" << P_min << ", " << P_max
                  << "] T=[" << T_min << ", " << T_max
                  << "] e_int=[" << e_int_min << ", " << e_int_max
                  << "] neg_rho=" << n_neg_rho << " neg_P=" << n_neg_P << std::endl;
        std::cout << "#   c[ENERGY]@ig=" << uCPU(ig_x, ig_y)[ENERGY]
                  << " rho*e_int=" << q1[DENSITY] * (q1[PRESSURE]/((GAMMA_DEFAULT-1.0)*(q1[DENSITY]+1e-30)))
                  << std::endl;
      }
    }

/*
std::stringstream filename1;
// 文件名固定为 "RESULT.txt"，不再动态生成
filename1 << solver.outputDirectory << "RESULT.txt";

std::ofstream outFile1;

// 第一次打开文件时以追加模式打开
outFile1.open(filename1.str().c_str(), std::ios::app);

outFile1.precision(8);
//outFile1 << "i" << " " << "j" << " " << "p" << ' ' << "rho" << ' ' << "t" << " " <<  "u" << " " <<  "v" << " " <<  "LA0" << " " <<  "LA1" << std::endl;
//for (int j = uCPU.activeNy()/5; j < uCPU.activeNy(); j += uCPU.activeNy()/5) {
  int j = ex_y * cell_y;  
  for (int i = 0; i < uCPU.activeNx(); ++i) {
      
      real q[NUMBER_VARIABLES];
      conservativeToPrimitive(uCPU(i, j), q);
      //outFile1 << std::fixed << uCPU(i, j)[DENSITY] << " " <<  uCPU(i, j)[YVELOCITY] / uCPU(i, j)[DENSITY] << std::endl;
      outFile1 << std::fixed << i << " " << j << " " << q[PRESSURE] << ' ' << q[DENSITY] << ' ' << q[PRESSURE]/q[DENSITY] << " " <<  q[XVELOCITY] << " " <<  q[YVELOCITY] << " " <<  q[SPECIES_INDEX(SP_CH4)] << " " <<  q[SPECIES_INDEX(SP_O2)] << std::endl;
    }
// 关闭文件流，确保数据写入完成
outFile1.close(); 
}
*/
      if (true) {
        std::vector<ImageOutputs> outputs;

        outputs.push_back((ImageOutputs){"pressure", PRESSUREFOO, HUE, 20000.0, 60000.0});
		//outputs.push_back((ImageOutputs){"sootfoil", PMAXPLOT, CUBEHELIX, 5.0, 50.0});
		//outputs.push_back((ImageOutputs){"sootfoil", PMAXPLOT, GREYSCALE, 10.0, 80.0});
		//outputs.push_back((ImageOutputs){"xvelocity", XVELOCITYPLOT, HUE, -1.0, 1.0});
    outputs.push_back((ImageOutputs){"velocity", VELOCITYPLOT, HUE, 0, 100.0});
        // outputs.push_back((ImageOutputs){"pressure2", PRESSUREFOO, HUE, 1.0, 50.0});
        outputs.push_back((ImageOutputs){"density", DENSITY, HUE, 0.05, 0.5});
       // outputs.push_back((ImageOutputs){"frac", FRACTIONPLOT, HUE, 0.0, 1.0});
	    //outputs.push_back((ImageOutputs){"LAMBDA0", YCH4PLOT, HUE, 0.0, 1.0});
        outputs.push_back((ImageOutputs){"YCH4", YCH4PLOT, HUE, 0.0, 0.142857});
        //outputs.push_back((ImageOutputs){"YO2", YO2PLOT, HUE, 0.0, 0.25});
        //outputs.push_back((ImageOutputs){"YH2O", YH2OPLOT, HUE, 0.0, 0.15});
        //outputs.push_back((ImageOutputs){"YCO2", YCO2PLOT, HUE, 0.0, 0.15});
        //outputs.push_back((ImageOutputs){"YOH", YOHPLOT, HUE, 0.0, 0.01});
        outputs.push_back((ImageOutputs){"schlieren", SCHLIEREN, GREYSCALE, 0.0, 1.0});
        outputs.push_back((ImageOutputs){"temperature", TEMPERATUREPLOT, HUE, 293, 800.0});

        for (std::vector<ImageOutputs>::iterator iter = outputs.begin(); iter != outputs.end(); ++iter) {
          std::stringstream filename;
          filename << solver.outputDirectory << (*iter).prefix << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".png";
          saveFrame(*solver.u, (*iter).plotVariable, (*iter).colourMode, filename.str().c_str(), (*iter).min, (*iter).max, & (*solver.fluxes)(0, 0, 0));
        }
      }
    }
  } while ((status = solver.step()) != Solver::FINISHED && !halt);

  times(&endTimes);
  clock_gettime(CLOCK_REALTIME, &endClock);
  const double wallTime = (endClock.tv_sec - startClock.tv_sec) + (endClock.tv_nsec - startClock.tv_nsec) * 1e-9;

  std::cout << "CPU time, wall= " << std::setprecision(2) << std::fixed << wallTime << "s, user=" << (endTimes.tms_utime - endTimes.tms_utime) << "s, sys=" << (endTimes.tms_stime - endTimes.tms_stime) << "s.  Time for {fluxes=" << solver.getTimeFluxes() * 1e-3 << "s, sources=" << solver.getTimeSourcing() * 1e-3 << "s, reduction=" << solver.getTimeReducing() * 1e-3 << "s, adding=" << solver.getTimeAdding() * 1e-3 << "s}" << std::endl;
  return 0;
}


