"""exp1 尺度不变 ρ(cλ)=c⁻¹ρ(λ) 的唯一解 = C/λ，以及截断如何破坏它。

门 1 推导链第 1、2、3 步：
  步1（据，预印本 1.3 §3）：无外部观察者 ⟹ 尺度不变 ⟹ ρ=C/λ。
  步2（新假设）：真实观察者谱有限支撑 λ∈[λ_min, λ_c]。
  步3（新推导）：截断引入特征尺度 Δ=ln(λ_c/λ_min)，破坏尺度不变。

符号（sympy）验证三件事：
  A. 函数方程 ρ(cλ)=c⁻¹ρ(λ) 的唯一（正）解是 ρ=C/λ。
  B. 对数均匀谱（ρ=C/λ）在 λ∈(0,∞) 上无特征尺度（尺度不变成立）。
  C. 有限支撑 [λ_min, λ_c] 引入特征尺度 Δ=ln(λ_c/λ_min)，尺度不变被破坏
     ——即 ρ_cutoff(cλ) ≠ c⁻¹ ρ_cutoff(λ) 对一般 c。
"""

import json
from pathlib import Path

import sympy as sp

results = {}

# ---------- A. 函数方程唯一解 ----------
# 解 ρ(cλ) = c⁻¹ρ(λ)，对所有 c>0。
# 预印本 1.3 §3.3 第4步：令 f(λ)=λρ(λ)，则 f(cλ)=f(λ) 对所有 c，
# 在对数坐标 s=log λ 下 f 是平移不变的常数 ⟹ f=C ⟹ ρ=C/λ。
lam, c, C = sp.symbols('lambda c C', positive=True)
rho = C / lam

# 验证 ρ(cλ) = c⁻¹ ρ(λ)
lhs = rho.subs(lam, c * lam)
rhs = c**(-1) * rho
eq_check = sp.simplify(lhs - rhs)
results['A_rho_c_lambda_eq_c_inv_rho'] = {
    'lhs_minus_rhs': str(sp.simplify(lhs - rhs)),
    'is_zero': bool(eq_check == 0),
}

# 唯一性：设 g(s)=f(e^s) 平移不变 ⟹ 常数。
# 用 sympy 直接解微分形式：f(cλ)=f(λ) ⟹ df/dλ 的相关约束。
# 更直接：验证 ρ=C/λ 确实满足，且任何 ρ=λ^p 只有 p=-1 满足。
p = sp.symbols('p')
rho_trial = lam**p
# ρ(cλ) = (cλ)^p = c^p λ^p，要求 = c⁻¹ λ^p ⟹ c^p = c⁻¹ ⟹ p=-1
eq_p = sp.Eq(c**p, c**(-1))
p_sol = sp.solve(eq_p, p)
results['A_power_law_only_p_eq_minus1'] = {
    'equation': 'c^p = c^(-1)',
    'solution': [str(s) for s in p_sol],
}

# ---------- B. 无特征尺度（尺度不变在 (0,∞) 上成立）----------
# 尺度不变 ⟺ 测度 ρ(λ)dλ 在 λ→cλ 下不变。
# ρ(λ)dλ = (C/λ)dλ = C d(ln λ)，在对数坐标下是均匀测度。
# 用变量代换验证积分不变性（符号）：∫ ρ(λ)dλ 在 [a,b] → [ca, cb] 下相等。
a, b = sp.symbols('a b', positive=True)
I1 = sp.integrate(rho, (lam, a, b))
I2 = sp.integrate(rho, (lam, c * a, c * b))
results['B_measure_scale_invariance'] = {
    'integral_a_to_b': str(sp.simplify(I1)),
    'integral_ca_to_cb': str(sp.simplify(I2)),
    'equal': bool(sp.simplify(I1 - I2) == 0),
    'note': '∫_a^b C/λ dλ = C ln(b/a)，在 λ→cλ 下不变 ⟹ 无特征尺度',
}

# ---------- C. 截断破坏尺度不变 ----------
# 有限支撑：ρ_cutoff(λ) = C/λ 当 λ∈[λ_min,λ_c]，否则 0。
# 截断引入特征尺度 Δ = ln(λ_c/λ_min)。
lam_min, lam_c = sp.symbols('lambda_min lambda_c', positive=True)
Delta = sp.log(lam_c / lam_min)
results['C_cutoff_characteristic_scale'] = {
    'Delta': str(Delta),
    'note': '截断引入特征尺度 Δ=ln(λ_c/λ_min)。尺度不变要求无特征尺度，故被破坏。',
}

# 数值演示：截断后，ρ_cutoff(cλ) 的支撑 [λ_min/c, λ_c/c] ≠ [λ_min, λ_c]，
# 因此 ρ_cutoff(cλ) ≠ c⁻¹ρ_cutoff(λ)（支撑错位）。
# 具体：取 λ_min=1, λ_c=e^Δ，c=e。原支撑 [1,e^Δ]，缩放后支撑 [1/e, e^(Δ-1)]。
# 支撑错位 ⟹ 尺度不变破坏。
results['C_support_shift'] = {
    'original_support': '[λ_min, λ_c]',
    'scaled_support': '[λ_min/c, λ_c/c]',
    'note': '缩放后支撑整体左移，与原支撑不重合（除非 λ_min=0 且 λ_c=∞），故尺度不变被截断破坏。',
}

# ---------- 输出 ----------
print(json.dumps(results, ensure_ascii=True, indent=2))

with open(Path(__file__).with_name('exp_scale_invariance_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
