"""exp4 门2 物理 Haldane 项（次近邻虚跃迁）：破 T 的物理扰动落在哪一类？

门 2 推导链第 4 步。

exp3 已坐实判据：破物理时间反演 T=K 的虚扰动 iM，Kramers 存活 ⟺ JM=−MJ。
本实验问：**真实的物理破 T 扰动（Haldane 型次近邻虚跃迁）落在哪一类？**

Haldane 质量（蜂窝格子标准）：次近邻跃迁带虚相位 t₂·e^{±iφ}，破时间反演。
在 π 磁通方格子上的对应：次近邻（对角）跃迁带虚相位 i·t₂·(符号)，
形成破 T 的复跳跃。

数值检查：
  A. 构造 π 磁通 D + Haldane 次近邻虚跃迁 H = D + i·t₂·M_haldane，
     看全偶重数（Kramers）是否存活。
  B. 计算 M_haldane 与 J 的对易/反对易关系，确定它落在 exp3 判据的哪一类。
"""

import json
from pathlib import Path

import numpy as np

results = {}


def toroidal_D(n_per_dim=4):
    n = n_per_dim ** 2
    D = np.zeros((n, n), float)
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            idx = n_per_dim * i + j
            jr = (j + 1) % n_per_dim
            D[idx, n_per_dim * i + jr] += 1.0
            D[n_per_dim * i + jr, idx] += 1.0
            id_ = (i + 1) % n_per_dim
            s = (-1.0) ** j
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
    n = D.shape[0]
    evals, evecs = np.linalg.eigh(D)
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
    return evecs @ J_block @ evecs.T


def all_even(H, tol):
    evals = np.linalg.eigvalsh(H)
    width = float(np.max(evals) - np.min(evals))
    scale = width if width > 1e-12 else 1.0
    mults = cluster_multiplicities(np.sort(evals), tol * scale)
    return bool(np.all(mults % 2 == 0)), mults


def haldane_M(n_per_dim):
    """次近邻（对角）虚跃迁：M 实反对称，H_haldane = i·t₂·M。

    方格子次近邻：(i,j)->(i+1,j+1) 与 (i,j)->(i+1,j-1)（对角）。
    破 T 的 Haldane 型：对角跃迁带虚相位，方向相关符号（形成环流）。
    M 需实反对称（iM 厄米）。
    """
    n = n_per_dim ** 2
    M = np.zeros((n, n), float)
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            idx = n_per_dim * i + j
            # 次近邻 1: (i,j) -> (i+1, j+1)，符号 +1
            i1 = (i + 1) % n_per_dim
            j1 = (j + 1) % n_per_dim
            n1 = n_per_dim * i1 + j1
            M[idx, n1] += 1.0
            M[n1, idx] -= 1.0  # 反对称
            # 次近邻 2: (i,j) -> (i+1, j-1)，符号 -1（破 T 的方向环流）
            j2 = (j - 1) % n_per_dim
            n2 = n_per_dim * i1 + j2
            M[idx, n2] -= 1.0
            M[n2, idx] += 1.0  # 反对称
    return M


n_per_dim = 4
n = n_per_dim ** 2
D = toroidal_D(n_per_dim)
J = build_real_J(D)

M_haldane = haldane_M(n_per_dim)
results['haldane_M_check'] = {
    'M_antisym_err': float(np.max(np.abs(M_haldane + M_haldane.T))),
    'M_real': bool(np.all(np.abs(M_haldane.imag) < 1e-12)) if np.iscomplexobj(M_haldane) else True,
}

# M_haldane 与 J 的关系：分解成反对易部分 + 对易部分
anti_part = (M_haldane + J @ M_haldane @ J) / 2  # 满足 JM=-MJ
comm_part = (M_haldane - J @ M_haldane @ J) / 2  # 满足 JM=+MJ
anti_norm = float(np.linalg.norm(anti_part))
comm_norm = float(np.linalg.norm(comm_part))
results['haldane_J_decomposition'] = {
    'anti_commuting_part_norm': anti_norm,
    'commuting_part_norm': comm_norm,
    'note': '若 anti>>comm，则 Haldane 主要落在「反对易」类，Kramers 应存活；反之破坏',
}

for t2 in [0.1, 0.3, 0.5, 1.0]:
    H = D + 1j * t2 * M_haldane
    even, mults = all_even(H, 1e-8)
    results[f't2_{t2}'] = {'all_even_mult': even, 'mults': mults.tolist()}

results['conclusion'] = (
    '看 Haldane 次近邻虚跃迁（真实的破 T 扰动）落在 exp3 判据的哪一类：'
    '若它主要与 J 反对易，则 spinless π 磁通 + Haldane 仍保持 Kramers 简并，'
    '这是框架独有的、可在真实材料里检验的预言。'
)

print(json.dumps(results, ensure_ascii=True, indent=2))
with open(Path(__file__).with_name('exp_kramers_haldane_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
