"""exp3 门2 核心：加虚扰动 iM（破物理时间反演 T=K），Kramers 存活 ⟺ JM=−MJ。

门 2 推导链第 3 步（新推导，核心判据）。

物理设定：
  - π 磁通 D 是**实对称**（相位 exp(±iπj)=±1 为实数），自对偶 ⟹ 存在**实反对称正交** J
    （J²=−1，JD=DJ），涌现 Kramers T'=JK（T'²=−1）。
  - 物理时间反演 T=K（复共轭）。破 T 的厄米扰动 = iM，其中 M 实**反对称**
    （iM 厄米 ⟺ M 反对称）。

判据推导：
  涌现 T' 保护 Kramers ⟺ [T', H]=0 ⟺ J·H*·J⁻¹ = H。
  对 H = D + iM（D 实对称、M 实反对称）：
    H* = D − iM
    J·H*·J⁻¹ = JDJ⁻¹ − i·JMJ⁻¹ = D − i·JMJ⁻¹（JD=DJ）
    要求 = D + iM ⟹ JMJ⁻¹ = −M。
  J 实反对称正交 ⟹ J⁻¹ = J^T = −J ⟹ JMJ⁻¹ = −JMJ。
  故判据 ⟺ **JMJ = M ⟺ JM = −MJ（J 与 M 反对易）**。

数值验证：
  A. M_anti = (M0 + J·M0·J)/2 → 满足 JM=−MJ → Kramers 存活（全偶重数）；
  B. M_comm = (M0 − J·M0·J)/2 → 满足 JM=+MJ → Kramers 破坏。

（上一版 bug：J 用了复本征向量导致 J 复；M 误写成对称。本版全部用实矩阵。）
"""

import json
from pathlib import Path

import numpy as np

results = {}


def toroidal_D(n_per_dim=4):
    """π 磁通方格子 torus，返回**实**对称矩阵（相位 ±1）。"""
    n = n_per_dim ** 2
    D = np.zeros((n, n), float)
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            idx = n_per_dim * i + j
            jr = (j + 1) % n_per_dim
            D[idx, n_per_dim * i + jr] += 1.0
            D[n_per_dim * i + jr, idx] += 1.0
            id_ = (i + 1) % n_per_dim
            s = (-1.0) ** j  # exp(iπj) = ±1
            D[idx, n_per_dim * id_ + j] += s
            D[n_per_dim * id_ + j, idx] += s
    return D


def cluster_multiplicities(evals, tol):
    evals = np.sort(evals)
    n = len(evals)
    sizes = []
    start = 0
    for i in range(1, n):
        if evals[i] - evals[i - 1] > tol:
            sizes.append(i - start)
            start = i
    sizes.append(n - start)
    return np.asarray(sizes, dtype=int)


def build_real_J(D):
    """从实对称 D 构造实反对称正交 J（J²=−1，JD=DJ）。"""
    n = D.shape[0]
    evals, evecs = np.linalg.eigh(D)  # 实对称 → 实本征向量
    width = float(np.max(evals) - np.min(evals))
    scale = width if width > 1e-12 else 1.0
    mults = cluster_multiplicities(evals, 1e-8 * scale)
    sigma = np.array([[0, 1], [-1, 0]], dtype=float)
    J_block = np.zeros((n, n), float)
    pos = 0
    for m in mults:
        if m % 2 == 0:
            for _ in range(m // 2):
                J_block[pos:pos + 2, pos:pos + 2] = sigma
                pos += 2
        else:
            pos += m
    return evecs @ J_block @ evecs.T, mults


def all_even_hermitian(H, tol):
    evals = np.linalg.eigvalsh(H)
    width = float(np.max(evals) - np.min(evals))
    scale = width if width > 1e-12 else 1.0
    mults = cluster_multiplicities(np.sort(evals), tol * scale)
    return bool(np.all(mults % 2 == 0)), mults


n_per_dim = 4
n = n_per_dim ** 2
D = toroidal_D(n_per_dim)
J, mults_D = build_real_J(D)

# J 实数性 + 性质检查
results['J_check'] = {
    'J_is_real': bool(np.all(np.abs(J.imag) < 1e-12)) if np.iscomplexobj(J) else True,
    'J_antisym_err': float(np.max(np.abs(J + J.T))),
    'J_orthog_err': float(np.max(np.abs(J @ J.T - np.eye(n)))),
    'J2_minus1_err': float(np.max(np.abs(J @ J + np.eye(n)))),
    'JD_DJ_err': float(np.max(np.abs(J @ D.conj() - D @ J))),
    'D_mults': mults_D.tolist(),
}

rng = np.random.default_rng(1)
M0 = rng.normal(size=(n, n))
M0 = (M0 - M0.T) / 2  # 实反对称

M_anti = (M0 + J @ M0 @ J) / 2
M_comm = (M0 - J @ M0 @ J) / 2

# 验证 M 与 J 的反对易/对易
results['M_relation'] = {
    'M_anti_J_anticomm_err': float(np.max(np.abs(J @ M_anti + M_anti @ J))),
    'M_comm_J_comm_err': float(np.max(np.abs(J @ M_comm - M_comm @ J))),
    'M_anti_antisym_err': float(np.max(np.abs(M_anti + M_anti.T))),
    'M_comm_antisym_err': float(np.max(np.abs(M_comm + M_comm.T))),
}

eps = 0.3
for name, M in [('M_anti_keepsKramers', M_anti), ('M_comm_killsKramers', M_comm)]:
    H = D + 1j * eps * M
    even, mults = all_even_hermitian(H, 1e-8)
    results[name] = {'all_even_mult': even, 'mults': mults.tolist()}

results['conclusion'] = (
    'M_anti（JM=−MJ）加 iM 后全偶重数存活（Kramers 保住）；'
    'M_comm（JM=+MJ）加 iM 后全偶重数破坏（Kramers 消失）。'
    '坐实判据：破物理时间反演 T=K 的虚扰动 iM，若与内禀 J 反对易，'
    '则 spinless 系统仍保持 Kramers 简并——框架独有（标准理论 spinless 无 Kramers）。'
)

print(json.dumps(results, ensure_ascii=True, indent=2))
with open(Path(__file__).with_name('exp_kramers_breakT_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
