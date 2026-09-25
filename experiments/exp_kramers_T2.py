"""exp1 门2 核心区分（符号）：涌现 T'=JK（T'²=−1）vs 物理时间反演 T=K（T²=+1）。

门 2 推导链第 1 步。

关键物理（据，预印本 1.1 定理 1 + exp7）：
  - spinless 费米子的物理时间反演是 T=K（复共轭），T²=+1（K²=1）。
    spinless 无自旋 ⟹ 无 Kramers 简并（Kramers 定理要求 T²=−1，需自旋 1/2）。
  - 但 π 磁通偶数尺寸「自对偶」⟹ 存在反幺 T'=JK（J 实反对称正交，J²=−1），
    T'² = JK·JK = J²K² = −1 ⟹ 涌现 Kramers 简并。
  - 这个 T'²=−1 来自 π 磁通（磁平移反对易 T_xT_y=−T_yT_x），不是外部时间反演。

符号（sympy）验证三件事：
  A. 反幺算子 T=UK（U 幺正，K 复共轭）的 T² 符号 = U·conj(U) 的符号。
  B. spinless 物理 T=K：U=I ⟹ T²=+1（无 Kramers）。
  C. 涌现 T'=JK：U=J（实反对称，J²=−1）⟹ T'²=−1（有 Kramers）。
  结论：spinless 系统「本不该有 Kramers」，涌现 T'²=−1 是 π 磁通给的、框架独有的来源。
"""

import json
from pathlib import Path

import sympy as sp

results = {}

# ---------- A. 反幺 T=UK 的 T² 符号 ----------
# T=UK，T² = UKUK = U·conj(U)·K² = U·conj(U)（因 K²=1）。
# 所以 T² 的符号 = U·conj(U) 的特征。
# U 幺正 ⟹ conj(U)=U^{-T}（转置逆）。T² = U·U^{-T}。
# 对实 U：conj(U)=U，T²=U²。
# 对反对称正交 U（U^T=-U, U²=-I）：T² = U·U = -I（因 U 实）。

# 用 2x2 具体矩阵验证
# 物理 T=K：U=I
I2 = sp.Matrix([[1, 0], [0, 1]])
T2_physical = I2 * I2  # U·conj(U)，U=I 实
results['A_physical_T'] = {
    'U': 'I（spinless，无自旋）',
    'T2': str(T2_physical),
    'sign': 'T²=+1 ⟹ 无 Kramers（Kramers 定理要求 T²=−1）',
}

# 涌现 T'=JK：U=J=[[0,1],[-1,0]]（实反对称正交，J²=−1）
J = sp.Matrix([[0, 1], [-1, 0]])
T2_emergent = J * J  # J 实 ⟹ conj(J)=J ⟹ T'² = J·J
results['B_emergent_Tprime'] = {
    'J': '[[0,1],[-1,0]]（实反对称正交，J²=-I）',
    'Tprime2': str(T2_emergent),
    'sign': 'T\'²=−I ⟹ 有 Kramers 简并',
}

# ---------- C. 验证 J 是唯一（至符号）实反对称 2x2 满足 J²=−I ----------
# （据，预印本 1.1 定理 5）
a = sp.symbols('a')
Jgen = sp.Matrix([[0, a], [-a, 0]])
J2 = Jgen * Jgen
sol = sp.solve(sp.Eq(J2, -sp.eye(2)), a)
results['C_J_uniqueness'] = {
    'general_antisym_2x2': '[[0,a],[-a,0]]',
    'J2': str(J2),
    'condition_J2_eq_minusI_solution_a': [str(s) for s in sol],
    'note': 'a=±1 是唯一解 ⟹ J=[[0,±1],[∓1,0]] 唯一（至符号），号差与 SU(2) 同源',
}

# ---------- D. Kramers 定理的核心：T² 决定简并结构 ----------
# T²=−1 ⟹ 本征值全偶重数（Kramers 对）；T²=+1 ⟹ 不强制。
# 这是「自对偶 ⟺ 全偶重数」的代数根（预印本 1.1 定理 1）。
results['D_kramers_essence'] = {
    'T2_minus1': '全偶重数（Kramers 对，ψ 与 Tψ 正交简并）',
    'T2_plus1': '不强制偶重数（无 Kramers 保护）',
    'note': 'spinless 物理 T²=+1 ⟹ 本无 Kramers；π 磁通涌现 T\'²=−1 ⟹ 强加 Kramers',
}

print(json.dumps(results, ensure_ascii=True, indent=2))
with open(Path(__file__).with_name('exp_kramers_T2_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
