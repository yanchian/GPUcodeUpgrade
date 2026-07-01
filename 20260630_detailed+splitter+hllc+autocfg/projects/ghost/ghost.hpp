/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
//#define GHOST
#define REACTIVE
#pragma once
#include "core.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <deque>
#include <queue>
#include <ctime>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <stdio.h>
#include <stddef.h>
#include <signal.h>
#include "boost/thread.hpp"
#include "Matrix.hpp"
#include "MatrixOperations.hpp"

#define USE_GL

// ============================================================
// 详细甲烷-空气反应机理 (C/H/O 物种 + N2)
// 物种数: 42 个 C/H/O 活性物种 + N2 惰性稀释剂 = 43
// ============================================================
#define NUM_SPECIES 43
#define DIMENSIONS 2
#define GC 2

#include "grid.hpp"

template<Processor P>
struct Mesh {
  typedef Mesh2D<real, P, NUMBER_VARIABLES> type;
};

template<Processor P>
struct LevelSet {
  typedef Mesh2D<real, P, 2> type;
};

StridedCell<real, NUMBER_VARIABLES> typedef Cell;

// ============================================================
// 物理常数与参考状态 (有量纲)
// ============================================================

// 通用气体常数 [J/(mol·K)]
const real R_UNIV = 8.314462618;

// 参考状态: 初始压力 20kPa, 初始温度 293K
const real P_REF = 20000.0;     // 20 kPa
const real T_REF = 293.0;       // 293 K

// 默认初始反应物比例 (可调)
// CH4 + 3 O2 (摩尔比 1:3, 当量比 phi = 2/3, 贫燃)
// 对应质量分数: CH4=16/(16+96)=0.142857, O2=96/(16+96)=0.857143
const real Y_CH4_INIT = 0.142857;  // 初始 CH4 质量分数
const real Y_O2_INIT  = 0.857143;  // 初始 O2 质量分数
const real Y_N2_INIT  = 0.0;       // 初始 N2 质量分数 (无 N2)
const real PHI_INIT   = 0.666667;  // 初始当量比

// Reynodls number
const real RE_NUM = 839.0489;

// Prandlt number 
const real PR_NUM = 1.0;

// Schmidt number 
const real SC_NUM = 1.0;

// 默认比热比 (用于 Riemann 求解器，空气近似)
const real GAMMA_DEFAULT = 1.4;

// 参考粘度 [Pa·s]
// mu0 用于扩散通量: μ = mu0 * (P/ρ)^s, 其中 P/ρ = R_mix * T
// CH4+3O2 混合物: R_mix≈296.9 J/(kg·K), μ_ref≈1.75e-5 Pa·s @ 293K
// mu0 = μ_ref / (R_mix*T_ref)^s = 1.75e-5 / (296.9*293)^0.76 = 1.75e-5 / 5633 = 3.11e-9
const real mu0 = 3.11e-9;

// 参考动力粘度 [Pa·s] (用于时间步长估计, CH4+3O2 @ 293K)
const real MU_REF = 1.75e-5;

// Sutherland 粘度温度指数 (空气 ≈ 0.76)
const real s = 0.76;

// ============================================================
// 有量纲长度参数
// 圆盘直径 = 90 mm, 计算域 = 100 mm × 100 mm
// 只需修改 N_CELLS_X/Y 和 L_DOMAIN_X/Y, cell_x/y 自动推导
// ============================================================
const int   N_CELLS_X = 2000;         // X方向网格数 (调试减半)
const int   N_CELLS_Y = 2000;         // Y方向网格数
const real L_DOMAIN_X = 0.1;         // 计算域X边长 [m] (100 mm)
const real L_DOMAIN_Y = 0.1;         // 计算域Y边长 [m] (100 mm)
const real D_DISK   = 0.09;         // 圆盘直径 [m] (90 mm)
const real R_DISK   = 0.045;        // 圆盘半径 [m] (45 mm)
const real cell_x   = (real)N_CELLS_X / L_DOMAIN_X;  // X方向每米网格点数 (自动推导)
const real cell_y   = (real)N_CELLS_Y / L_DOMAIN_Y;  // Y方向每米网格点数 (自动推导)

///////////////////////////////////////////////////////////////////////////////////////////*/
const real PI = 3.14159265359;
const real Tol = 1.0e-6;//tolerance
const real Plot_step = 0.5;//numbers of steps between outputs. Playing a similar role as interval in .cfg file

real Counter = 1.0;
real OutputCounter = 6.0;
real Check = 1.0;
const real Check3 = 1.0;

//Inflow parameter
const real m = 0; //Mach number
const real theta = 0; //turn angle
const real wedge_x = 0; //starting point of the wedge
//
//Triangle geometry
const real A_x = 0;
const real B_x = 0;
const real AB_y = L_DOMAIN_Y;
const real theta_C = 45;//in degrees

//circle boundary (有量纲)
const real l_bump = L_DOMAIN_X * 0.5 * 1.414;//Length on the wedge surface (center to wedge tip)
const real r_bump = 0.045;  // 圆盘半径 [m] (直径 = 90 mm)

//explosion center (有量纲)
const real ex_x = L_DOMAIN_X * 0.5;
const real ex_y = L_DOMAIN_Y * 0.5;
const real R_IGNITION = 0.0027;   // 点火区半径 [m] (2.7 mm, 等比缩放: 0.015*0.045/0.2484)
const real ex_r = R_IGNITION;     // 点火区半径
// 爆燃初始条件: 点火区压力与环境压力同量级, 避免激波
const real ex_p = 12;              // 点火区压力倍数: P_ignition = ex_p * P_REF = 40 kPa
const real ex_rho = 0.3;         // 点火区密度 [kg/m^3], T ≈ 2712 K (理想气体)

// ============================================================
// 详细化学反应机理数据 (device-side constant memory)
// 由 parse_therm.ps1 从 THERM.DAT 自动生成
// ============================================================

// 分子量 [kg/mol]
__device__ __constant__ real MOL_WT[NUM_SPECIES];

// NASA 多项式温度切换点 [K]
__device__ __constant__ real NASA_TMID[NUM_SPECIES];

// NASA 多项式系数: [species][14]
// 布局: low(a1..a7), high(a1..a7)
__device__ __constant__ real NASA_COEFFS[NUM_SPECIES][14];

// ============================================================
// 热力学函数 (device)
// ============================================================

// 单物种 Cp, H, S 计算 (NASA 7-coefficient polynomial)
__device__ __forceinline__ void nasaCpHRTS(int isp, real T, real& cp, real& h, real& s) {
    int idx = (T < NASA_TMID[isp]) ? 0 : 7;
    const real* a = &NASA_COEFFS[isp][idx];
    
    real cpR = a[0];
    real hRT = a[0];
    real Tk = T;
    for (int k = 1; k < 5; k++) {
        cpR += a[k] * Tk;
        hRT += a[k] * Tk / (real)(k + 1);
        Tk *= T;
    }
    hRT += a[5] / T;
    
    real logT = log(T);
    real sR = a[0] * logT;
    Tk = T;
    for (int k = 1; k < 5; k++) {
        sR += a[k] * Tk / (real)k;
        Tk *= T;
    }
    sR += a[6];
    
    cp = cpR * R_UNIV;   // J/(mol·K)
    h = hRT * R_UNIV * T; // J/mol
    s = sR * R_UNIV;      // J/(mol·K)
}

// 混合气体内能 + 比热容 + 焓 [J/kg], [J/(kg·K)], [J/kg] — 融合计算, 省两次物种循环
__device__ __forceinline__ void mixtureThermoECH(const real* Y, real T, real& e, real& cp, real& h) {
    e = 0.0;
    cp = 0.0;
    h = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
        if (Y[i] > 0.0) {
            real cp_i, h_i, s_i;
            nasaCpHRTS(i, T, cp_i, h_i, s_i);
            real invM = 1.0 / MOL_WT[i];
            real Yi_invM = Y[i] * invM;
            cp += Yi_invM * cp_i;
            real h_i_invM = h_i * invM;
            real R_i = R_UNIV * invM;
            e += Y[i] * (h_i_invM - R_i * T);
            h += Yi_invM * h_i;
        }
    }
}

// 混合气体比热容 [J/(kg·K)]
__device__ __forceinline__ real mixtureCp(const real* Y, real T) {
    real cp = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
        if (Y[i] > 0.0) {
            real cp_i, h_i, s_i;
            nasaCpHRTS(i, T, cp_i, h_i, s_i);
            cp += Y[i] * cp_i / MOL_WT[i];  // cp_i J/(mol·K), MOL_WT kg/mol
        }
    }
    return cp;
}

// 混合气体焓 [J/kg]
__device__ __forceinline__ real mixtureEnthalpy(const real* Y, real T) {
    real h = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
        if (Y[i] > 0.0) {
            real cp_i, h_i, s_i;
            nasaCpHRTS(i, T, cp_i, h_i, s_i);
            h += Y[i] * h_i / MOL_WT[i];
        }
    }
    return h;
}

// 混合气体生成焓 [J/kg] (at T_ref = 298.15 K)
// 用于化学热释放能量更新: Q = h_form_old - h_form_new
__device__ __forceinline__ real mixtureFormationEnthalpy(const real* Y) {
    const real T_ref = 298.15;
    real h = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
        if (Y[i] > 0.0) {
            real cp_i, h_i, s_i;
            nasaCpHRTS(i, T_ref, cp_i, h_i, s_i);
            h += Y[i] * h_i / MOL_WT[i];
        }
    }
    return h;
}

// 混合气体内能 [J/kg]
__device__ __forceinline__ real mixtureInternalEnergy(const real* Y, real T) {
    real e = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
        if (Y[i] > 0.0) {
            real cp_i, h_i, s_i;
            nasaCpHRTS(i, T, cp_i, h_i, s_i);
            real R_i = R_UNIV / MOL_WT[i];  // J/(kg·K)
            e += Y[i] * (h_i / MOL_WT[i] - R_i * T);
        }
    }
    return e;
}

// 混合气体平均分子量 [kg/mol]
__device__ __forceinline__ real mixtureMolWt(const real* Y) {
    real mwInv = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
        if (Y[i] > 0.0) {
            mwInv += Y[i] / MOL_WT[i];
        }
    }
    return 1.0 / (mwInv + 1e-30);
}

// 混合气体常数 [J/(kg·K)]
__device__ __forceinline__ real mixtureRgas(const real* Y) {
    return R_UNIV / mixtureMolWt(Y);
}

// 从内能求解温度 (Newton 迭代) — 使用融合 e+cp+h 计算, 返回焓供能量更新
__device__ __forceinline__ real solveTemperature(const real* Y, real e, real T_guess, real& h_mix) {
    real T = T_guess;
    real R_mix = mixtureRgas(Y);  // 混合气体常数在迭代中不变, 提到外面
    for (int iter = 0; iter < 10; iter++) {
        real e_mix, cp_mix;
        mixtureThermoECH(Y, T, e_mix, cp_mix, h_mix);
        real f = e_mix - e;
        real cv = cp_mix - R_mix;
        real dT = f / (cv + 1e-30);
        T -= dT;
        if (fabs(dT / (T + 1e-30)) < 1e-6) break;
    }
    return T;
}

// 从 Cell 提取质量分数
__device__ __forceinline__ void getMassFractions(const Cell u, real* Y) {
    real sumY = 0.0;
    for (int k = 0; k < NUM_SPECIES - 1; k++) {
        Y[k] = u[SPECIES_INDEX(k)];
        sumY += Y[k];
    }
    Y[NUM_SPECIES - 1] = 1.0 - sumY;
}

// 从 Cell 计算温度
__device__ __forceinline__ real temperature(const Cell u) {
    real rho = u[DENSITY];
    real e_kin = 0.5 * (u[XVELOCITY] * u[XVELOCITY] + u[YVELOCITY] * u[YVELOCITY]);
    real e_int = (u[ENERGY] / (rho + 1e-30)) - e_kin;
    
    real Y[NUM_SPECIES];
    getMassFractions(u, Y);
    
    real h_dummy;
    return solveTemperature(Y, e_int, T_REF, h_dummy);
}

// 从 Cell 计算压力
__device__ __forceinline__ real pressureFromCell(const Cell u) {
    real T = temperature(u);
    real rho = u[DENSITY];
    real Y[NUM_SPECIES];
    getMassFractions(u, Y);
    real R_mix = mixtureRgas(Y);
    return rho * R_mix * T;
}

// 计算混合比热比 gamma = cp/cv
__device__ __forceinline__ real gammaMix(const Cell u) {
    real Y[NUM_SPECIES];
    getMassFractions(u, Y);
    real T = temperature(u);
    real cp = mixtureCp(Y, T);
    real R_mix = mixtureRgas(Y);
    real cv = cp - R_mix;
    return cp / (cv + 1e-30);
}



////////////////////////////////////////////////////////////////////////////////
//const int len_inert = 500; //in terms of cells (multiplied by resolution)
real density_new = 1.0;
const real Delta_rho = 0.0;
const real frame_interface = 1.0;

const real Output_step = 0.5;
real Plot_counter = 10.0;
const int Skip_lines = 5;

/////////////////////////////////////////////////////////////////
__device__ __host__ __forceinline__ real mu() {
  return MU_REF;
}

__device__ __host__ __forceinline__ real gamma(const Cell u) {
  return GAMMA_DEFAULT;
}

__device__ __host__ __forceinline__ real gamma() {
  return GAMMA_DEFAULT;
}

__device__ __host__ __forceinline__ real p0(const Cell u) {
  return 0.0;
}
__device__ __host__ __forceinline__ real p0() {
  return 0.0;
}

// 包含机理数据 (host-side)
#include "mechanism_data.ipp"

#include "boundaryconditions.hpp"
#include "flux.hpp"
#include "diffusionflux.hpp"
#include "wavespeed.hpp"
#include "find_minimum.hpp"
#include "HLLC.hpp"
#include "HLLE.hpp"
#include "Solver.hpp"
#include "initialconditions.hpp"
#include "render.hpp"
#include "opengl.hpp"
#ifdef GHOST
#include "ghostfluid.hpp"
#endif
#include "source.hpp"

#ifdef REACTIVE
#include "source.hpp"
#include "shockdetect.hpp"
#endif

struct ImageOutputs {
  std::string prefix;
  int plotVariable;
  ColourMode colourMode;
  real min;
  real max;
};

#include "SDF/BoundaryMesh.hpp"
#include "SDF/Polyhedron.cu"
#include "SDF/ConnectedEdge.cu"
#include "SDF/Edge.cu"
#include "SDF/Face.cu"
#include "SDF/Vertex.cu"
#include "SDF/ConnectedFace.cu"
#include "SDF/ConnectedVertex.cu"
#include "SDF/ScanConvertiblePolygon.cu"
#include "SDF/ScanConvertiblePolyhedron.cu"
#include "SDF/BoundaryMesh.cu"

#include "kernels/boundaryconditions.ipp"
#include "kernels/flux.ipp"
#include "kernels/diffusionflux.ipp"
#ifdef GHOST
#include "kernels/ghostfluid.ipp"
#endif

#ifdef REACTIVE
//#include "kernels/source.ipp"
#include "kernels/shockdetect.ipp"
#endif

#include "kernels/source.ipp"
#include "kernels/HLLC.ipp"
#include "kernels/HLLE.ipp"
#include "kernels/wavespeed.ipp"
#include "kernels/find_minimum.ipp"
