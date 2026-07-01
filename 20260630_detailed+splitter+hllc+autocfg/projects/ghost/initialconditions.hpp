/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
  Modified: Dimensional initial conditions for detailed CH4-O2 mechanism
  Reference: P0 = 20 kPa, T0 = 293 K
  Mixture: CH4 + 3 O2 (molar ratio 1:3, phi = 2/3)
*/
#pragma once

// 物种索引宏 (与 mechanism_data.ipp 中 SPECIES_NAMES_HOST 的顺序一致)
#define SP_H     0
#define SP_H2    1
#define SP_O     2
#define SP_O2    3
#define SP_OH    4
#define SP_HO2   5
#define SP_H2O   6
#define SP_H2O2  7
#define SP_CO    8
#define SP_CO2   9
#define SP_HCO  10
#define SP_CH3  11
#define SP_CH4  12
#define SP_CH2O 13
#define SP_N2   42

// 点火区自由基初始质量分数 (用于加速起爆)
#define Y_H_INIT   0.001   // H 原子
#define Y_O_INIT   0.001   // O 原子
#define Y_OH_INIT  0.003   // OH 自由基
#define Y_RAD_TOTAL (Y_H_INIT + Y_O_INIT + Y_OH_INIT)  // 自由基总质量分数 = 0.005

// 点火区偏移量 (等比缩放: 0.0319 * 0.045/0.2484 = 0.00578 m)
#define IGNITION_OFFSET 0.0058

__global__ void setInitialConditions(Mesh<GPU>::type u) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x,
            j = blockIdx.y * blockDim.y + threadIdx.y;

  if (u.active(i, j)) {
    const real x = u.x(i), y = u.y(j);

    real rho, P, v_x = 0.0, v_y = 0.0;

    // 单点点火: 中心高温区 (爆燃模式, 压力与环境同量级避免激波)
    // 双点点火已注释, 调试用单点
    // bool inIgnition = (pow(x - ex_x - IGNITION_OFFSET, 2) + pow(y - ex_y, 2) <= pow(ex_r, 2) || 
    //                    pow(x - ex_x + IGNITION_OFFSET, 2) + pow(y - ex_y, 2) <= pow(ex_r, 2));
    bool inIgnition = (pow(x - ex_x, 2) + pow(y - ex_y, 2) <= pow(ex_r, 2));
    
    if (inIgnition) {
      P   = ex_p * P_REF;   // 40 kPa (2× 环境压力), 避免强激波
      rho = ex_rho;         // 0.05 kg/m³, T ≈ P/(rho*R_mix) ≈ 2712 K
    }
    else {
      P = P_REF;  // 20 kPa reference pressure

      // 混合物平均分子量: MW_mix = 1 / sum(Y_k / MW_k)  [kg/mol]
      // CH4+3O2: Y_CH4=0.142857, Y_O2=0.857143
      real invMW = Y_CH4_INIT / MOL_WT[SP_CH4] + 
                   Y_O2_INIT  / MOL_WT[SP_O2];
      real MW_mix = 1.0 / (invMW + 1e-30);

      // 理想气体: rho = P * MW_mix / (R_univ * T)
      rho = P_REF * MW_mix / (R_UNIV * T_REF);
    }

    // 总能量: E = P/(gamma-1) + 0.5*rho*u^2
    // 仅包含显热（不含生成焓），由 Riemann 求解器（γ=1.4 常数 EOS）使用
    // 化学焓变化通过 source.ipp 中的 c[ENERGY] += rho*dh 补充
    real energy = P / (GAMMA_DEFAULT - 1.0) + 0.5 * rho * (v_x*v_x + v_y*v_y);

    // 存储初始守恒变量
    u(i, j, DENSITY) = rho;
    u(i, j, XMOMENTUM) = v_x * rho;
    u(i, j, YMOMENTUM) = v_y * rho;
    u(i, j, ENERGY) = energy;

    // 设置初始物种质量分数 (使用上面已定义的 inIgnition)
    // bool inIgnition = (pow(x - ex_x - IGNITION_OFFSET, 2) + pow(y - ex_y, 2) <= pow(ex_r, 2) || 
    //                    pow(x - ex_x + IGNITION_OFFSET, 2) + pow(y - ex_y, 2) <= pow(ex_r, 2));

    for (int k = 0; k < NUM_SPECIES - 1; k++) {
      real Yk = 0.0;

      if (inIgnition) {
        // 点火区: CH4 + O2 + 微量自由基 (H, O, OH)
        // 自由基总计 0.5%, CH4和O2等比缩减
        real scale = 1.0 - Y_RAD_TOTAL;  // 0.995
        if (k == SP_CH4) {
          Yk = Y_CH4_INIT * scale;
        } else if (k == SP_O2) {
          Yk = Y_O2_INIT * scale;
        } else if (k == SP_H) {
          Yk = Y_H_INIT;
        } else if (k == SP_O) {
          Yk = Y_O_INIT;
        } else if (k == SP_OH) {
          Yk = Y_OH_INIT;
        }
      } else {
        // 非点火区: CH4 + O2, 无自由基
        if (k == SP_CH4) {
          Yk = Y_CH4_INIT;
        } else if (k == SP_O2) {
          Yk = Y_O2_INIT;
        }
      }

      u(i, j, SPECIES_INDEX(k)) = Yk;
    }

    u(i, j, ISSHOCK) = 0.0;
    u(i, j, PMAX) = P;
  }
}