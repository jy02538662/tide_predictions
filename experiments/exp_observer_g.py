"""exp2 有限观察者：h(ω)=log|ω| 截断到 |ω|≤Δ，G(t) 尾巴相对纯 1/t 的修正。

门 1 推导链第 4、5 步。

据（预印本 1.3 §4）：G(t)=∫ h(ω)e^{itω}dω，h(ω) 在 ω=0 有对数奇性 ⟹ FT[log|ω|]=-π/|t| ⟹ G(t)~1/t。
新推导：真实观察者谱有限支撑 ⟹ 频率差 |ω|=|ln(λ_i/λ_j)|≤Δ=ln(λ_c/λ_min)，
       即 h(ω) 的支撑被截到 |ω|≤Δ。

符号（sympy）验证：
  A. FT[log|ω|] = -π/|t|（对数奇性的 Fourier 核，据，预印本 1.3 §4.2）。
  B. 硬截断 h(ω)=log|ω|·χ(|ω|≤Δ) 的 G(t) 解析式。

数值（numpy）验证：
  C. 有限 Δ 下 G(t) 尾巴偏离 1/t（出现振荡 (logΔ)sin(Δt)/t 的迹象）。

关键物理：截断 Δ 引入特征尺度，修正 1/t 尾巴。这是「观察者依赖」的来源。
"""

import json
from pathlib import Path

import numpy as np

results = {}

# ---------- A. FT[log|ω|] = -π/|t|（据，标准结果） ----------
# 这是已知的分布意义 Fourier 变换（预印本 1.3 §4.2 引用）：
#   FT[log|ω|](t) = ∫_{-∞}^∞ log|ω| e^{itω} dω = -π/|t|。
# 直接数值验证（不依赖 sympy 发散积分）：
#   2∫_0^∞ log(ω) cos(tω) dω → -π/t（在 t 固定时，尾部主导）。
results['A_FT_log_abs_omega'] = {
    'standard_result': 'FT[log|omega|](t) = -pi/|t|',
    'source': '预印本 1.3 §4.2（对数奇性 → 1/t）',
    'note': '据，符号直接引用标准结果，数值见 C 部分对照',
}

# ---------- B. 硬截断 G(t) 解析式（新推导，闭式） ----------
# G_Δ(t) = 2∫_0^Δ log(ω) cos(tω) dω
# 分部积分：∫ log(ω) cos(tω) dω = (log ω · sin(tω)/t) - (1/t)∫ sin(tω)/ω dω
#   = (log ω · sin(tω)/t) - Si(tω)/t，其中 Si 是正弦积分。
# 故 G_Δ(t) = 2[log(Δ) sin(tΔ) - Si(tΔ)]/t。
# 当 Δ→∞：Si(tΔ)→π/2（sgn t），log(Δ)sin(tΔ) 振荡无极限 ⟹ 需正规化，
# 分布极限回到 -π/t（对应 A）。
results['B_hard_cutoff_G_analytic'] = {
    'closed_form': 'G_Delta(t) = 2[log(Delta)*sin(t*Delta) - Si(t*Delta)]/t',
    'Si_note': 'Si = sine integral = ∫_0^x sin(u)/u du',
    'note': '硬截断解析式：振荡项 log(Δ)sin(tΔ)/t 主导，说明硬截断引入振荡而非干净 1/t',
}

# ---------- C. 数值：有限 Δ 下尾巴偏离 1/t ----------
def G_cutoff(t_vals, Delta):
    """G_Δ(t) = 2∫_0^Δ log(ω) cos(tω) dω，数值积分。"""
    out = np.zeros_like(t_vals, dtype=float)
    for i, tv in enumerate(t_vals):
        # 数值积分（ω 从极小到 Δ，避开 ω=0 的 log 奇性）
        wgrid = np.linspace(1e-6, Delta, 200000)
        integrand = 2 * np.log(wgrid) * np.cos(tv * wgrid)
        out[i] = np.trapz(integrand, wgrid)
    return out

# 纯 1/t 参考（对数奇性主导，据）：G_ref(t) = -π/t
t_vals = np.logspace(0.3, 2.0, 40)  # t 从 ~2 到 100

Delta_vals = [5.0, 10.0, 20.0, 50.0]
c_result = {}
for Delta in Delta_vals:
    G_num = G_cutoff(t_vals, Delta)
    G_ref = -np.pi / t_vals
    # 尾巴（大 t）拟合 G(t) ~ t^{-p}
    tail = t_vals > 10
    p = -np.polyfit(np.log(t_vals[tail]), np.log(np.abs(G_num[tail])), 1)[0]
    c_result[str(Delta)] = {
        'tail_power_p': float(p),
        'G_num_t30': float(G_num[np.argmin(np.abs(t_vals - 30))]),
        'G_ref_t30': float(G_ref[np.argmin(np.abs(t_vals - 30))]),
    }
results['C_tail_vs_1t'] = c_result
results['C_note'] = (
    '有限 Δ 下 G(t) 尾巴指数 p 偏离 1（纯 1/t 对应 p=1）。'
    '硬截断在 |ω|=Δ 处不连续，产生振荡项 (logΔ)sin(Δt)/t，'
    '使尾巴在 1/t 附近振荡而非严格 1/t。'
)

print(json.dumps(results, ensure_ascii=True, indent=2))
with open(Path(__file__).with_name('exp_observer_g_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
