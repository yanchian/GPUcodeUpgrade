/*
  Detailed CH4-Air reaction mechanism (C/H/O species only)
  Based on A. Konnov's mechanism v0.4, filtered to C/H/O species
  Inert N2 is retained as bath gas
*/

#pragma once

// ============================================================
// 物种定义 - 共 40 个 C/H/O 物种 + N2 (惰性)
// ============================================================
#define NUM_SPECIES 41

enum SpeciesIndex {
  // H2/O2 基础体系
  S_H, S_H2, S_O, S_O2, S_OH, S_HO2, S_H2O, S_H2O2,
  // C1 物种
  S_CO, S_CO2, S_HCO, S_CH3, S_CH4, S_CH2O, S_CH2, S_CH3O,
  S_CH2OH, S_CH, S_C, S_SCH2, S_CH3O2, S_CH3O2H,
  // C2 物种
  S_C2H6, S_C2H5, S_C2H4, S_C2H3, S_C2H2, S_C2H, S_CH2CO, S_HCCO,
  S_CH3OH, S_CH3HCO, S_CH3CO, S_CH2HCO, S_C2H4O, S_C2H5O, S_C2H5O2,
  S_C2H5O2H, S_C2, S_C2O, S_CH3CO3, S_CH3CO3H,
  // N2 (惰性稀释剂)
  S_N2
};

// 分子量 (g/mol)
__device__ __constant__ real MOL_WT[NUM_SPECIES] = {
  1.00797,    // H
  2.01594,    // H2
  15.99940,   // O
  31.99880,   // O2
  17.00737,   // OH
  33.00677,   // HO2
  18.01534,   // H2O
  34.01474,   // H2O2
  28.01055,   // CO
  44.00995,   // CO2
  29.01852,   // HCO
  15.03506,   // CH3
  16.04303,   // CH4
  30.02649,   // CH2O
  14.02709,   // CH2
  31.03446,   // CH3O
  31.03446,   // CH2OH
  13.01912,   // CH
  12.01115,   // C
  14.02709,   // SCH2
  47.03386,   // CH3O2
  48.04183,   // CH3O2H
  30.07012,   // C2H6
  29.06215,   // C2H5
  28.05418,   // C2H4
  27.04621,   // C2H3
  26.03824,   // C2H2
  25.03027,   // C2H
  42.03764,   // CH2CO
  41.02967,   // HCCO
  32.04243,   // CH3OH
  44.05358,   // CH3HCO
  43.04561,   // CH3CO
  43.04561,   // CH2HCO
  44.05358,   // C2H4O
  45.06155,   // C2H5O
  61.06095,   // C2H5O2
  62.06892,   // C2H5O2H
  24.02230,   // C2
  40.02170,   // C2O
  75.04441,   // CH3CO3
  76.05238,   // CH3CO3H
  28.01340    // N2
};

// ============================================================
// NASA 多项式系数 (7-coefficient, 200-6000K)
// Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
// H/RT = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
// S/R  = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
// 每个物种有 2 组 (低温 + 高温), 共 14 个系数
// ============================================================
#define NASA_COEFFS_PER_SPECIES 14
#define NASA_LOW_RANGE  0
#define NASA_HIGH_RANGE 7

// 通用气体常数 [J/(mol*K)]
#define RU 8.314462618

// 温度切换点 (K)
__device__ __constant__ real NASA_TMID[NUM_SPECIES];

// NASA 多项式系数: [species][14]
// 布局: low(a1..a7), high(a1..a7)
__device__ __constant__ real NASA_COEFFS[NUM_SPECIES][NASA_COEFFS_PER_SPECIES];

// ============================================================
// 反应结构定义
// ============================================================

// 反应类型
enum ReactionType {
  ARRHENIUS,        // 标准 Arrhenius: k = A * T^b * exp(-Ea/RT)
  THIRD_BODY,       // 三体碰撞: k = [M] * A * T^b * exp(-Ea/RT)
  LINDEMANN,        // Lindemann 形式
  TROE,             // Troe 形式
  DUPLICATE          // 重复反应 (与其他反应共用)
};

// 单个反应的数据
struct ReactionData {
  real A;              // 指前因子 [cm^3,mol,s]
  real b;              // 温度指数
  real Ea;             // 活化能 [cal/mol]
  ReactionType type;   // 反应类型
  short nReactants;    // 反应物数量
  short nProducts;     // 产物数量
  short reactants[3];  // 反应物物种索引 (最多 3 个)
  short products[4];   // 产物物种索引 (最多 4 个)
  short stoichReactants[3];  // 反应物化学计量系数
  short stoichProducts[4];   // 产物化学计量系数
  // 三体效率 (仅 THIRD_BODY, LINDEMANN, TROE)
  bool hasThirdBodyEff;
  short thirdBodySpecies[8];   // 增强三体效率的物种
  real thirdBodyEff[8];        // 对应效率
  // 低压极限 (LINDEMANN, TROE)
  real lowA, lowB, lowEa;
  // Troe 参数
  real troeAlpha, troeT3, troeT1, troeT2;
};

#define MAX_REACTIONS 350
__device__ __constant__ ReactionData REACTIONS[MAX_REACTIONS];
__device__ __constant__ int NUM_REACTIONS;

// ============================================================
// 热力学函数 (device)
// ============================================================
__device__ __forceinline__ void nasaCpHRTS(
    int isp, real T, real& cp, real& h, real& s) {
  
  int idx = (T < NASA_TMID[isp]) ? NASA_LOW_RANGE : NASA_HIGH_RANGE;
  const real* a = &NASA_COEFFS[isp][idx];

  real Ti = 1.0;
  real cpR = a[0];
  real hRT = a[0];
  for (int k = 1; k < 5; k++) {
    Ti *= T;
    cpR += a[k] * Ti;
    hRT += a[k] * Ti / (k + 1);
  }
  hRT += a[5] / T;

  real sR = a[0] * log(T);
  for (int k = 1; k < 5; k++) {
    sR += a[k] * pow(T, k) / k;
  }
  sR += a[6];

  cp = cpR * RU;
  h = hRT * RU * T;
  s = sR * RU;
}

// 混合气体比热容 [J/(kg*K)]
__device__ __forceinline__ real mixtureCp(const real* Y, real T) {
  real cp = 0.0;
  for (int i = 0; i < NUM_SPECIES; i++) {
    if (Y[i] > 0.0) {
      real cp_i, h_i, s_i;
      nasaCpHRTS(i, T, cp_i, h_i, s_i);
      cp += Y[i] * cp_i / (MOL_WT[i] * 1e-3);  // cp_i in J/(mol*K), convert to J/(kg*K)
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
      h += Y[i] * h_i / (MOL_WT[i] * 1e-3);
    }
  }
  return h;
}

// 混合气体平均分子量 [kg/mol]
__device__ __forceinline__ real mixtureMolWt(const real* Y) {
  real mwInv = 0.0;
  for (int i = 0; i < NUM_SPECIES; i++) {
    if (Y[i] > 0.0) {
      mwInv += Y[i] / (MOL_WT[i] * 1e-3);
    }
  }
  return 1.0 / (mwInv + 1e-30);
}

// 混合气体比热比
__device__ __forceinline__ real mixtureGamma(const real* Y, real T) {
  real cp = mixtureCp(Y, T);
  real mw = mixtureMolWt(Y);
  real cv = cp - RU / mw;
  return cp / (cv + 1e-30);
}

// ============================================================
// 反应速率计算
// ============================================================

// 计算单个 Arrhenius 反应速率 [cm^3/mol/s] 或 [1/s]
__device__ __forceinline__ real arrheniusRate(const ReactionData& rxn, real T) {
  real k = rxn.A * pow(T, rxn.b) * exp(-rxn.Ea / (1.9872037 * T));  // R = 1.987 cal/(mol*K)
  return k;
}

// 计算三体浓度 [mol/cm^3]
__device__ __forceinline__ real thirdBodyConcentration(
    const real* conc, const ReactionData& rxn) {
  real sum = 0.0;
  if (rxn.hasThirdBodyEff) {
    for (int i = 0; i < 8 && rxn.thirdBodySpecies[i] >= 0; i++) {
      sum += conc[rxn.thirdBodySpecies[i]] * (rxn.thirdBodyEff[i] - 1.0);
    }
  }
  // 所有物种默认效率为 1.0
  for (int i = 0; i < NUM_SPECIES; i++) {
    sum += conc[i];
  }
  return sum;
}

// 计算反应进度 (单个反应)
__device__ __forceinline__ real reactionProgress(
    const ReactionData& rxn, real T, const real* conc) {
  
  real kf = arrheniusRate(rxn, T);

  // 计算反应物浓度乘积
  real reactProd = 1.0;
  for (int i = 0; i < rxn.nReactants; i++) {
    int sp = rxn.reactants[i];
    if (sp >= 0 && sp < NUM_SPECIES) {
      real c = conc[sp];
      real exponent = rxn.stoichReactants[i];
      if (exponent > 1) {
        reactProd *= pow(c, exponent);
      } else {
        reactProd *= c;
      }
    }
  }

  real rate = kf * reactProd;

  // 三体/压力修正
  if (rxn.type == THIRD_BODY || rxn.type == LINDEMANN || rxn.type == TROE) {
    real M = thirdBodyConcentration(conc, rxn);
    rate *= M;
  }

  // Lindemann/Troe 低压修正
  if (rxn.type == LINDEMANN || rxn.type == TROE) {
    real k0 = rxn.lowA * pow(T, rxn.lowB) * exp(-rxn.lowEa / (1.9872037 * T));
    real M = thirdBodyConcentration(conc, rxn);
    real Pr = k0 * M / (kf + 1e-60);
    
    if (rxn.type == TROE) {
      real logPr = log10(Pr + 1e-60);
      real c = -0.4 - 0.67 * log10(rxn.troeAlpha);
      real n = 0.75 - 1.27 * log10(rxn.troeAlpha);
      real d = 0.14;
      real Fcent = (1.0 - rxn.troeAlpha) * exp(-T / rxn.troeT3) 
                 + rxn.troeAlpha * exp(-T / rxn.troeT1) 
                 + exp(-rxn.troeT2 / T);
      real logF = log10(Fcent) / (1.0 + pow((logPr + c) / (n - d * (logPr + c)), 2));
      real F = pow(10.0, logF);
      real Finf = Pr / (1.0 + Pr);
      rate *= F * Finf;
    } else {
      // Lindemann
      real Finf = Pr / (1.0 + Pr);
      rate *= Finf;
    }
  }

  return rate;
}