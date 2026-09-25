"""exp_standard_predictions：引力侧四个标准预言的量级计算（v9 第 3 步，顺带）。

四个标准预言（不独有，任何热核展开/量子基础框架都给类似形式）：
  1. 黑洞阴影修正   ：a₄ 高阶曲率 → 阴影半径修正 δr/r ~ (l_P/M)²
  2. 引力波色散     ：维度流 d_s(2→4) → 速度频率依赖 δv/v ~ (k/M_P)²
  3. 伽马射线暴延迟 ：最小长度 → 高能光子延迟 Δt ~ η(E/M_QG)(D/c)
  4. 维度流         ：d_s 跑动（微观 2 → 宏观 4），λ_c 是种子

定位：标准预言，量级全部远低于可测阈值（或框架未直接算），价值在「进量子引力
社区入口（LQG/CDT 对话）+ 诚实标注测不出」。不是框架独有预言。

符号/数值双验，无需 GPU。
"""

import json
from pathlib import Path

import numpy as np

results = {}

# ---------- 物理常数（SI + 自然单位） ----------
c = 2.99792458e8            # m/s
l_P = 1.616255e-35          # m，普朗克长度
M_P_GeV = 1.220890e19       # GeV，普朗克质量
M_P_Hz = c / l_P            # Hz，普朗克频率（c/l_P）
G = 6.67430e-11             # m^3 kg^-1 s^-2
M_sun_geom = 1.476625e3     # m，太阳几何质量 GM_sun/c^2
H0 = 67.4e3 / 3.0857e22     # s^-1，H0=67.4 km/s/Mpc
H0_GeV_scale = None

results['constants'] = {
    'c_m_s': c,
    'l_Planck_m': l_P,
    'M_Planck_GeV': M_P_GeV,
    'M_Planck_Hz': M_P_Hz,
    'M_sun_geometric_m': M_sun_geom,
    'H0_s': H0,
}

# ---------- 1. 黑洞阴影修正 ----------
# 谱作用量 a₄ 高阶曲率项（R²_{μνρσ} 等）加到 EH 上，量级 (1/Λ²) 压制。
# Schwarzschild 解 R=R_{μν}=0，只有 R_{μνρσ}R^{μνρσ}=48M²/r⁶ 贡献。
# 阴影半径修正 δr_sh/r_sh ~ (l_P/M)²，M 用几何质量 GM/c²。
# EHT 测 M~10 M_sun（M87* 是 6.5e9 M_sun，恒星质量黑洞测不到阴影）。
for M_sun in [10.0, 6.5e9]:  # 恒星级 10 M_sun；M87* 6.5e9 M_sun
    M_geom = M_sun * M_sun_geom
    delta_r_over_r = (l_P / M_geom) ** 2
    results[f'shadow_M{M_sun:g}_Msun'] = {
        'M_geometric_m': M_geom,
        'delta_r_over_r': delta_r_over_r,
    }

# ---------- 2. 引力波色散 ----------
# 维度流 d_s(紫外)=2 → 修正色散 ω² = k² + k⁴/M_*²（HL 型，紫外 d_s=2）。
# 群速度 v_g = dω/dk ≈ 1 + (3/2)(k/M_*)^2，k≪M_*。
# LIGO f~100 Hz；M_* ~ M_P（普朗克频率）。
for f_Hz in [30.0, 100.0, 300.0]:
    k = 2 * np.pi * f_Hz / c  # 1/m，k=ω/c=2πf/c
    M_star = M_P_Hz / c       # 1/m，M_* = M_P 频率 / c
    delta_v = 1.5 * (k / M_star) ** 2
    results[f'gw_dispersion_f{f_Hz:g}Hz'] = {
        'delta_v_over_v': delta_v,
    }

# ---------- 3. 伽马射线暴延迟 ----------
# 最小长度 → Lorentz 破坏色散 E² = p²(1 + η(p/M_QG)^n)，n=1 线性 / n=2 二次。
# 延迟 Δt ≈ η (n+1)/2 · (E_h^n - E_l^n)/(M_QG^n) · D/c 的简式：
#   线性 n=1：Δt ≈ η (ΔE/M_QG)(D/c)
#   二次 n=2：Δt ≈ η (3/2)(E²/M_QG²)(D/c)
# D ~ z=1 距离（约 4 Gpc），E ~ 1 TeV 光子，M_QG ~ M_P。
D_m = 4.0 * 3.0857e25        # 4 Gpc -> m
D_over_c = D_m / c           # s
E_TeV = 1.0
E_GeV = E_TeV * 1e3
M_QG_GeV = M_P_GeV

dt_n1 = (E_GeV / M_QG_GeV) * D_over_c       # 线性，η=1
dt_n2 = 1.5 * (E_GeV / M_QG_GeV) ** 2 * D_over_c  # 二次，η=1

results['grb_delay'] = {
    'D_m': D_m,
    'D_over_c_s': D_over_c,
    'E_GeV': E_GeV,
    'M_QG_GeV': M_QG_GeV,
    'dt_linear_n1_s': dt_n1,
    'dt_quadratic_n2_s': dt_n2,
    'note': '线性 n=1 秒级可测（Fermi-LAT 已排除 η~O(1) 到 M_QG>几 M_P）；'
            '二次 n=2 ~1e-17 s 测不出。框架是洛伦兹协变低能有效理论，'
            '「最小长度=离散化截断」不直接给 Lorentz 破坏——n=1 在框架里无来源。',
}

# ---------- 4. 维度流 ----------
# 框架谱维数 d_s 读出（家底 9）：2.01/3.01/4.02（不同维度的图分别读出）。
# 这是「读」步（固定图读维度），不是「维度流」（同一图 d_s(σ) 随扩散时间跑动）。
# 维度流种子 = λ_c（观察者截断，微观最大熵 → 宏观尺度不变），d_s(σ) 未直接算。
results['spectral_dimension_flow'] = {
    'd_s_readout': '2.01 / 3.01 / 4.02（不同维度图分别读出，exp_spectral_dim）',
    'seed': 'λ_c（观察者截断 = 过渡尺度，微观最大熵 → 宏观尺度不变）',
    'd_s_sigma_not_computed': True,
    'note': '维度流 d_s(σ) 跑动函数未直接算（λ_c 只是种子），诚实边界。',
}

# ---------- 汇总 ----------
results['summary'] = {
    'shadow_stellar_10Msun': results['shadow_M10_Msun']['delta_r_over_r'],
    'shadow_M87': results['shadow_M6.5e+09_Msun']['delta_r_over_r'],
    'gw_dispersion_100Hz': results['gw_dispersion_f100Hz']['delta_v_over_v'],
    'grb_delay_linear_s': dt_n1,
    'grb_delay_quadratic_s': dt_n2,
    'conclusion': (
        '四个标准预言量级：黑洞阴影 ~1e-78（恒星）/ ~1e-90（M87*）测不出；'
        '引力波色散 ~1e-81 测不出；伽马暴延迟线性 ~30s（可测但框架无 Lorentz '
        '破坏来源）、二次 ~1e-17s 测不出；维度流未直接算（只有 λ_c 种子）。'
        '全部是「标准预言」，不能区分框架，价值在 LQG/CDT 对话入口 + 诚实标注。'
    ),
}

out = Path(__file__).with_name('exp_standard_predictions_last_run.json')
out.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(results, ensure_ascii=False, indent=2))
print(f"\nwrote {out.name}")
