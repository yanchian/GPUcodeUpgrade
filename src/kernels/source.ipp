/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
  Modified: Detailed CH4-Air reaction mechanism with point implicit integration
  Reference pressure 20kPa, reference temperature 293K
*/
#include "source.hpp"

// ============================================================
// 精简核心反应机理 (C/H/O 物种, 约 20 个关键反应)
// Arrhenius: k = A * T^b * exp(-Ea/RT)
// A: [cm^3,mol,s], Ea: [cal/mol], R: 1.987 [cal/(mol*K)]
// ============================================================

#define NREACTIONS 22
#define RU_CAL 1.987  // universal gas constant in cal/(mol*K)

// 反应物/产物物种索引 (在 mechanism_data.ipp 中定义的顺序)
// H=0, H2=1, O=2, O2=3, OH=4, HO2=5, H2O=6, H2O2=7, CO=8, CO2=9,
// HCO=10, CH3=11, CH4=12, CH2O=13, CH2=14, CH3O=15, CH2OH=16, CH=17, C=18,
// SCH2=19, CH3O2=20, CH3O2H=21, C2H6=22, C2H5=23, C2H4=24, ...

__device__ __forceinline__ void computeReactionRates(
    real T, const real* C,     // C: species concentrations [mol/cm^3]
    real* kf,                   // forward rate constants
    real* wdot                  // species net production rates [mol/(cm^3·s)]
) {
  for (int i = 0; i < NUM_SPECIES; i++) wdot[i] = 0.0;

  real Ti = 1.0 / T;
  real logT = log(T);

  // 预计算总浓度 M (用于第三体反应, 避免重复遍历)
  real M_total = 0.0;
  for (int i = 0; i < NUM_SPECIES; i++) M_total += C[i];

  // Reaction 1: H + O2 = O + OH
  {
    real A = 2.644e16, b = -0.6707, Ea = 17041.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[0] * C[3];  // S_H * S_O2
    wdot[0] -= rate;  wdot[3] -= rate;
    wdot[2] += rate;  wdot[4] += rate;
  }

  // Reaction 2: O + H2 = H + OH
  {
    real A = 4.589e4, b = 2.700, Ea = 6260.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[2] * C[1];  // S_O * S_H2
    wdot[2] -= rate;  wdot[1] -= rate;
    wdot[0] += rate;  wdot[4] += rate;
  }

  // Reaction 3: OH + H2 = H + H2O
  {
    real A = 1.734e8, b = 1.510, Ea = 3430.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[4] * C[1];  // S_OH * S_H2
    wdot[4] -= rate;  wdot[1] -= rate;
    wdot[0] += rate;  wdot[6] += rate;
  }

  // Reaction 4: OH + OH = O + H2O
  {
    real A = 3.973e4, b = 2.400, Ea = -2110.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[4] * C[4];  // S_OH * S_OH
    wdot[4] -= 2.0 * rate;
    wdot[2] += rate;  wdot[6] += rate;
  }

  // Reaction 5: H + O2 + M = HO2 + M
  {
    real A = 4.651e12, b = 0.440, Ea = 0.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real M = M_total;
    real rate = k * C[0] * C[3] * M;
    wdot[0] -= rate;  wdot[3] -= rate;
    wdot[5] += rate;  // S_HO2
  }

  // Reaction 6: H + HO2 = O + H2O
  {
    real A = 3.970e12, b = 0.000, Ea = 671.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[0] * C[5];  // S_H * S_HO2
    wdot[0] -= rate;  wdot[5] -= rate;
    wdot[2] += rate;  wdot[6] += rate;
  }

  // Reaction 7: H + HO2 = H2 + O2
  {
    real A = 6.620e13, b = 0.000, Ea = 2130.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[0] * C[5];  // S_H * S_HO2
    wdot[0] -= rate;  wdot[5] -= rate;
    wdot[1] += rate;  wdot[3] += rate;
  }

  // Reaction 8: OH + HO2 = H2O + O2
  {
    real A = 2.130e15, b = -0.576, Ea = 0.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[4] * C[5];  // S_OH * S_HO2
    wdot[4] -= rate;  wdot[5] -= rate;
    wdot[6] += rate;  wdot[3] += rate;
  }

  // Reaction 9: CO + OH = CO2 + H
  {
    real A = 4.760e7, b = 1.228, Ea = 70.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[8] * C[4];  // S_CO * S_OH
    wdot[8] -= rate;  wdot[4] -= rate;
    wdot[9] += rate;  wdot[0] += rate;
  }

  // Reaction 10: CH4 + H = CH3 + H2
  {
    real A = 5.470e7, b = 1.970, Ea = 11210.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[12] * C[0];  // S_CH4 * S_H
    wdot[12] -= rate;  wdot[0] -= rate;
    wdot[11] += rate;  wdot[1] += rate;
  }

  // Reaction 11: CH4 + OH = CH3 + H2O
  {
    real A = 5.720e6, b = 1.960, Ea = 2639.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[12] * C[4];  // S_CH4 * S_OH
    wdot[12] -= rate;  wdot[4] -= rate;
    wdot[11] += rate;  wdot[6] += rate;
  }

  // Reaction 12: CH3 + O = CH2O + H
  {
    real A = 5.540e13, b = 0.050, Ea = -136.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[11] * C[2];  // S_CH3 * S_O
    wdot[11] -= rate;  wdot[2] -= rate;
    wdot[13] += rate;  wdot[0] += rate;
  }

  // Reaction 13: CH2O + H = HCO + H2
  {
    real A = 2.300e10, b = 1.050, Ea = 3275.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[13] * C[0];  // S_CH2O * S_H
    wdot[13] -= rate;  wdot[0] -= rate;
    wdot[10] += rate;  wdot[1] += rate;
  }

  // Reaction 14: CH2O + OH = HCO + H2O
  {
    real A = 3.430e9, b = 1.180, Ea = -447.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[13] * C[4];  // S_CH2O * S_OH
    wdot[13] -= rate;  wdot[4] -= rate;
    wdot[10] += rate;  wdot[6] += rate;
  }

  // Reaction 15: HCO + M = H + CO + M
  {
    real A = 1.870e17, b = -1.000, Ea = 17000.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real M = M_total;
    real rate = k * C[10] * M;  // S_HCO * M
    wdot[10] -= rate;
    wdot[0] += rate;  wdot[8] += rate;  // S_CO
  }

  // Reaction 16: HCO + O2 = CO + HO2
  {
    real A = 1.345e13, b = 0.000, Ea = 400.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[10] * C[3];  // S_HCO * S_O2
    wdot[10] -= rate;  wdot[3] -= rate;
    wdot[8] += rate;  wdot[5] += rate;
  }

  // Reaction 17: CH3 + H(+M) = CH4(+M)  [recombination]
  {
    real A = 1.270e16, b = -0.630, Ea = 383.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real M = M_total;
    real rate = k * C[11] * C[0] * M;  // S_CH3 * S_H * M
    wdot[11] -= rate;  wdot[0] -= rate;
    wdot[12] += rate;  // S_CH4
  }

  // Reaction 18: H + OH + M = H2O + M
  {
    real A = 4.400e22, b = -2.000, Ea = 0.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real M = M_total;
    real rate = k * C[0] * C[4] * M;
    wdot[0] -= rate;  wdot[4] -= rate;
    wdot[6] += rate;
  }

  // Reaction 19: H + H + M = H2 + M
  {
    real A = 1.780e18, b = -1.000, Ea = 0.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real M = M_total;
    real rate = k * C[0] * C[0] * M;
    wdot[0] -= 2.0 * rate;
    wdot[1] += rate;
  }

  // Reaction 20: O + O + M = O2 + M
  {
    real A = 6.165e15, b = -0.500, Ea = 0.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real M = M_total;
    real rate = k * C[2] * C[2] * M;
    wdot[2] -= 2.0 * rate;
    wdot[3] += rate;
  }

  // Reaction 21: H2O2 + OH = H2O + HO2
  {
    real A = 1.000e12, b = 0.000, Ea = 0.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[7] * C[4];  // S_H2O2 * S_OH
    wdot[7] -= rate;  wdot[4] -= rate;
    wdot[6] += rate;  wdot[5] += rate;
  }

  // Reaction 22: CH4 + O2 = CH3 + HO2  (初始化反应, 产生第一批自由基)
  // GRI-Mech 3.0: A=4.0e13, b=0.0, Ea=56900 cal/mol
  {
    real A = 4.000e13, b = 0.000, Ea = 56900.0;
    real k = A * exp(b * logT - Ea * Ti / RU_CAL);
    real rate = k * C[12] * C[3];  // S_CH4 * S_O2
    wdot[12] -= rate;  wdot[3] -= rate;
    wdot[11] += rate;  wdot[5] += rate;  // S_CH3, S_HO2
  }
}

// ============================================================
// 点隐式源项积分
// dY_k/dt = wdot_k * M_k / rho  (mass fraction rate)
// 使用对角近似: Y_k^{n+1} = (Y_k^n + dt * P_k) / (1 + dt * D_k)
// ============================================================
template<int blockDimX, int blockDimY>
__global__ void sources(Mesh<GPU>::type u, const real dt) {
  const int dimX = blockDimX,
            dimY = blockDimY;
  const int i = dimX * blockIdx.x + threadIdx.x;
  const int j = dimY * blockIdx.y + threadIdx.y;

  if (u.active(i, j)) {
    real c[NUMBER_VARIABLES];
    real p[NUMBER_VARIABLES];
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      c[k] = u(i, j, k);
    }

    conservativeToPrimitive(c, p);

    // 提取质量分数
    real Y[NUM_SPECIES];
    real sumY = 0.0;
    for (int k = 0; k < NUM_SPECIES - 1; k++) {
      Y[k] = p[SPECIES_INDEX(k)];
      sumY += Y[k];
    }
    Y[NUM_SPECIES - 1] = 1.0 - sumY;

    // 计算温度: 使用理想气体定律 T = P/(ρ*R_mix) = (γ-1)*e_int/R_mix
    // 与 Riemann 求解器 (γ=1.4 常数 EOS) 保持一致
    // 不使用 solveTemperature (含生成焓的热力学与显热 EOS 不兼容)
    real rho = p[DENSITY];
    real e_kin = 0.5 * (p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY]);
    real e_int = (c[ENERGY] / (rho + 1e-30)) - e_kin;
    real R_mix = mixtureRgas(Y);
    real T = (GAMMA_DEFAULT - 1.0) * e_int / (R_mix + 1e-30);

    // 只在高温区计算化学反应 (T > 1000K 或存在自由基)
    // 先检查当前值是否有效
    if (isnan(T) || isnan(rho) || rho < 1e-10) {
      // 状态无效, 跳过反应
      for (int k = 0; k < NUMBER_VARIABLES; k++) {
        u(i, j, k) = c[k];
      }
      return;
    }
    bool hasRadicals = (Y[0] > 1e-10 || Y[4] > 1e-10 || Y[2] > 1e-10);
    if (T > 1000.0 || hasRadicals) {

      // 计算物种浓度 [mol/cm^3]
      // rho: kg/m^3 -> g/cm^3 (×1e-3)
      // MOL_WT[k]: kg/mol -> g/mol (×1e3)
      // C = Y*rho_gcm3/(MOL_WT*1e3) = Y*rho*1e-6/MOL_WT [mol/cm^3]
      real C[NUM_SPECIES];
      real rho_gcm3 = rho * 1e-3;  // kg/m^3 -> g/cm^3
      for (int k = 0; k < NUM_SPECIES; k++) {
        if (Y[k] > 1e-30) {
          C[k] = Y[k] * rho_gcm3 / (MOL_WT[k] * 1e3);  // mol/cm^3
        } else {
          C[k] = 0.0;
        }
      }

      // 点隐式更新质量分数
      real dt_sub = dt;
      int nsub = 1;
      // 化学时间尺度: 自由基反应 ~1e-9s, 使用更多子步提高稳定性
      real tau_chem = 1e-9;
      if (dt > tau_chem) {
        nsub = max(1, (int)(dt / tau_chem));
        dt_sub = dt / nsub;
      }
      // 限制最大子步数为 500, 防止超慢
      if (nsub > 500) {
        nsub = 500;
        dt_sub = dt / nsub;
      }

      // 保存反应前生成焓, 用于反应后能量更新
      real h_form_old = mixtureFormationEnthalpy(Y);

      for (int isub = 0; isub < nsub; isub++) {
        // 重新计算反应速率 (因为组分已更新)
        real wdot_sub[NUM_SPECIES];
        computeReactionRates(T, C, nullptr, wdot_sub);

        real dY[NUM_SPECIES];
        for (int k = 0; k < NUM_SPECIES; k++) {
          real Mk = MOL_WT[k] * 1e3;  // kg/mol -> g/mol
          real omega_k = wdot_sub[k] * Mk / rho_gcm3;  // dY_k/dt [1/s]
          
          // NaN 检查
          if (isnan(omega_k)) omega_k = 0.0;
          
          // 点隐式: 分离产生项和消耗项
          real P_k = 0.0, D_k = 0.0;
          if (omega_k > 0.0) {
            P_k = omega_k;
          } else {
            D_k = -omega_k / (Y[k] + 1e-30);
          }
          
          dY[k] = (Y[k] + dt_sub * P_k) / (1.0 + dt_sub * D_k) - Y[k];
          if (isnan(dY[k])) dY[k] = 0.0;
        }

        // 更新组分
        for (int k = 0; k < NUM_SPECIES; k++) {
          Y[k] += dY[k];
          if (Y[k] < 0.0 || isnan(Y[k])) Y[k] = 0.0;
        }

        // 归一化
        real sumYnew = 0.0;
        for (int k = 0; k < NUM_SPECIES; k++) sumYnew += Y[k];
        if (sumYnew < 1e-30) {
          // 全部物种耗尽, 回退
          break;
        }
        real invSum = 1.0 / sumYnew;
        for (int k = 0; k < NUM_SPECIES; k++) Y[k] *= invSum;

        // 更新浓度
        for (int k = 0; k < NUM_SPECIES; k++) {
          if (Y[k] > 1e-30) {
            C[k] = Y[k] * rho_gcm3 / (MOL_WT[k] * 1e3);
          } else {
            C[k] = 0.0;
          }
        }

        // 子步内能量更新: 化学热释放 Q = h_form_old - h_form_new (正值为放热)
        R_mix = mixtureRgas(Y);
        real h_form_new = mixtureFormationEnthalpy(Y);
        real Q_chem = h_form_old - h_form_new;
        
        // 能量释放率限制: 每子步最多增加 2% 内能
        real max_Q = 0.02 * (e_int + 1e-10);
        if (Q_chem > max_Q) Q_chem = max_Q;
        if (Q_chem < -max_Q) Q_chem = -max_Q;
        if (isnan(Q_chem)) Q_chem = 0.0;
        
        c[ENERGY] += rho * Q_chem;
        e_int = c[ENERGY] / (rho + 1e-30) - e_kin;
        T = (GAMMA_DEFAULT - 1.0) * e_int / (R_mix + 1e-30);
        
        // 温度上限 5000K + NaN 保护
        const real T_MAX = 5000.0;
        if (T > T_MAX || isnan(T)) {
          T = T_MAX;
          e_int = T_MAX * R_mix / (GAMMA_DEFAULT - 1.0);
          c[ENERGY] = rho * (e_int + e_kin);
        }
        if (e_int < 0.0 || isnan(e_int)) {
          e_int = 1.0;
          c[ENERGY] = rho * (e_int + e_kin);
          T = (GAMMA_DEFAULT - 1.0) * e_int / (R_mix + 1e-30);
        }
        h_form_old = h_form_new;
      }
    }

    // 将质量分数写回
    for (int k = 0; k < NUM_SPECIES - 1; k++) {
      c[SPECIES_INDEX(k)] = Y[k];
    }

    // NaN 最终检查: 写回前确保所有守恒量有效
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (isnan(c[k]) || isinf(c[k])) {
        c[k] = 0.0;
      }
    }
    if (c[DENSITY] < 1e-10) c[DENSITY] = 1e-10;
    if (c[ENERGY] < 1e-10) c[ENERGY] = 1e-10;

    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      u(i, j, k) = c[k];
    }
  }
}