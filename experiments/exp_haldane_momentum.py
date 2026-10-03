"""
验证 Haldane M 的动量空间形式 + 为什么保 Kramers 简并（修复 FFT，用解析 + 数值直接验证）。

核心：caveat 里「M(k)=M(-k)」是错的。正确结论：
  1. M̂(k) = 4i cos(kx) sin(ky)（含 sin，奇函数 M̂(-k) = -M̂(k)）。
  2. iM̂(k) = -4 cos(kx) sin(ky)（实，厄米项）。
  3. 关键机制：iM̂ 在「谷间平移」k→k+(π,π) 下不变（iM̂(k+(π,π)) = iM̂(k)），
     而涌现 T' 的 Kramers 对正是谷间 (k, k+(π,π))，所以 Haldane 项对两者给相同修正，简并存活。
     不是 M(k)=M(-k)，是谷间平移不变。
"""
import numpy as np
import json
from pathlib import Path

results = {}

# ---- 解析推导（直接用解析，避开 FFT 约定）----
# M 平移不变，非零元（反对称）：M(1,1)=+1, M(-1,-1)=-1, M(1,-1)=-1, M(-1,1)=+1
# M̂(k) = Σ M(δ) e^{i k·δ}
#       = e^{i(kx+ky)} - e^{-i(kx+ky)} - e^{i(kx-ky)} + e^{-i(kx-ky)}
#       = 2i sin(kx+ky) - 2i sin(kx-ky) = 2i[2 cos(kx) sin(ky)] = 4i cos(kx) sin(ky)
def M_hat(kx, ky):
    return 4j * np.cos(kx) * np.sin(ky)

def iM_hat(kx, ky):  # iM̂(k)（实，厄米项的动量空间形式）
    return -4.0 * np.cos(kx) * np.sin(ky)

# ---- 1. 解析形式 ----
results['M_hat_form'] = {
    'M̂(k)': '4i cos(kx) sin(ky)',
    'iM̂(k)（厄米项）': '-4 cos(kx) sin(ky)（实）',
    '含 sin（奇函数）': True,
}

# ---- 2. 数值验证：随机 k 点，M̂ 解析 vs 直接算 M 的傅里叶 ----
n_per_dim = 6
n = n_per_dim ** 2
# 直接构造 M 实空间，再算傅里叶（用平移不变性，M̂(k) = Σ M(0->δ) e^{ikδ}）
# 从 (0,0) 出发的 δ 位移：M(1,1)=+1, M(1,-1)=-1（及反对称，但反对称项由 e^{ikδ}+e^{-ikδ} 自动含）
# 精确：M̂(k) = Σ_{δ 有序} M(δ) e^{ikδ}，M(δ) 反对称
def M_hat_numeric(kx, ky):
    s = 0j
    # 四个非零 δ（含反对称）
    for dx, dy, val in [(1,1,1),(-1,-1,-1),(1,-1,-1),(-1,1,1)]:
        s += val * np.exp(1j * (kx * dx + ky * dy))
    return s

rng = np.random.default_rng(0)
max_err = 0.0
for _ in range(100):
    kx = rng.uniform(-np.pi, np.pi)
    ky = rng.uniform(-np.pi, np.pi)
    a = M_hat(kx, ky)
    b = M_hat_numeric(kx, ky)
    max_err = max(max_err, abs(a - b))
results['M_hat_numeric'] = {
    '100 随机 k 点，M̂ 解析 vs 直接傅里叶 最大误差': f'{max_err:.2e}',
    'M̂(k) = 4i cos(kx) sin(ky) 坐实': bool(max_err < 1e-12),
}

# ---- 3. 奇函数：M̂(-k) = -M̂(k) ----
kx, ky = 0.7, 1.3
results['M_hat_odd'] = {
    'M̂(k)': f'{M_hat(kx,ky):.6f}',
    'M̂(-k)': f'{M_hat(-kx,-ky):.6f}',
    'M̂(-k) = -M̂(k)（奇函数）': f'{abs(M_hat(-kx,-ky) + M_hat(kx,ky)):.2e}',
    '结论': 'M̂ 是奇函数 M̂(-k)=-M̂(k)，不是偶函数 M(k)=M(-k)——用户对',
}

# ---- 4. 谷间平移不变：iM̂(k+(π,π)) = iM̂(k)（保简并的正确机制）----
kx, ky = 0.5, 0.8
val_k = iM_hat(kx, ky)
val_valley = iM_hat(kx + np.pi, ky + np.pi)
results['valley_invariance'] = {
    'iM̂(k)': f'{val_k:.6f}',
    'iM̂(k+(π,π))': f'{val_valley:.6f}',
    'iM̂(k+(π,π)) = iM̂(k)（谷间平移不变）': f'{abs(val_valley - val_k):.2e}',
    '解析验证': 'cos(kx+π)=-cos(kx), sin(ky+π)=-sin(ky) ⟹ iM̂(k+(π,π)) = -4(-cos kx)(-sin ky) = iM̂(k)',
    '结论': 'Haldane 项对谷间 Kramers 对 (k, k+(π,π)) 给相同修正，简并存活。正确机制是「谷间平移不变」，不是「M(k)=M(-k)」',
}

# ---- 5. 但注意：iM̂ 是奇函数（k→-k），所以 Haldane 确实破 k→-k 对称 ----
results['note'] = {
    'iM̂(-k) = -iM̂(k)': '奇函数，E(k)≠E(-k)，Haldane 破 k→-k 对称',
    '但 Kramers 对不是 (k,-k)': '涌现 T\'=JK 是谷间 k→k+(π,π)（把 Dirac 点 (±π/2,±π/2) 互映），不是 k→-k',
    '所以保护简并的是': 'iM̂ 的谷间平移不变性（iM̂(k+(π,π))=iM̂(k)），不是 k→-k 的偶性',
}

print(json.dumps(results, ensure_ascii=False, indent=2))
with open(Path(__file__).with_name('exp_haldane_momentum_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
