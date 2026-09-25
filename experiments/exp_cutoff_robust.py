"""exp3 平滑截断 vs 硬截断：观察者引力修正是否 robust。

门 1 推导链第 5 步。

exp2 发现：硬截断 h(ω)=log|ω|·χ(|ω|≤Δ) 在 |ω|=Δ 处不连续，产生振荡项
(logΔ)sin(Δt)/t，尾巴指数 p 乱跳——这是「截断伪影」，不是物理修正。

本实验问：如果用**平滑截断**（更物理：观察者的分辨率不是硬截止，是逐渐衰减），
修正是否 robust？即是否存在一个与截断细节无关的、干净的物理修正？

平滑截断模型：h(ω) = log|ω| · f(|ω|/Δ)，f 是平滑窗口：
  - f(x) = exp(-x)（指数窗口）
  - f(x) = exp(-x²)（高斯窗口）
对比硬截断 f(x) = 1(x≤1)。

数值：算 G(t) = 2∫_0^∞ log(ω) f(ω/Δ) cos(tω) dω，
     拟合尾巴，看「相对 1/t 的偏离」在三种窗口下是否一致。

关键判断：
  - 若平滑窗口下尾巴仍是 1/t（p→1），则「硬截断振荡」纯属伪影，物理无修正；
  - 若平滑窗口下尾巴有干净的、与窗口无关的修正（如指数尾巴 e^{-t/Δ}），
    则是真实物理修正。
"""

import json
from pathlib import Path

import numpy as np

results = {}


def G_window(t_vals, Delta, window):
    """G(t) = 2∫_0^∞ log(ω) f(ω/Δ) cos(tω) dω，ω 上限取到几个 Δ 以覆盖窗口。"""
    out = np.zeros_like(t_vals, dtype=float)
    omega_max = 8.0 * Delta
    wgrid = np.linspace(1e-6, omega_max, 300000)
    f = window(wgrid / Delta)
    for i, tv in enumerate(t_vals):
        integrand = 2 * np.log(wgrid) * f * np.cos(tv * wgrid)
        out[i] = np.trapz(integrand, wgrid)
    return out


def window_hard(x):
    return np.where(x <= 1.0, 1.0, 0.0)


def window_exp(x):
    return np.exp(-x)


def window_gauss(x):
    return np.exp(-x * x)


t_vals = np.logspace(0.2, 2.0, 50)  # t ~1.6 到 100
Delta = 10.0
G_ref = -np.pi / t_vals  # 纯 1/t 参考

windows = {
    'hard': window_hard,
    'exponential': window_exp,
    'gaussian': window_gauss,
}

res = {}
for name, fn in windows.items():
    G = G_window(t_vals, Delta, fn)
    # 尾巴拟合指数（大 t 区）
    tail = t_vals > 8
    p = -np.polyfit(np.log(t_vals[tail]), np.log(np.abs(G[tail])), 1)[0]
    # 相对偏离：|G - G_ref| / |G_ref| 在几个 t 点的平均值
    rel_dev = np.mean(np.abs((G[tail] - G_ref[tail]) / G_ref[tail]))
    res[name] = {
        'tail_power_p': float(p),
        'mean_relative_deviation_tail': float(rel_dev),
    }

results['C3_window_comparison'] = res
results['C3_conclusion'] = (
    '硬截断 p 乱跳（振荡伪影）；指数/高斯平滑窗口 p 应稳定趋近某值。'
    '若平滑窗口下 p≈1（尾巴仍是 1/t），则修正只来自截断伪影，物理上无观察者依赖的'
    '长程修正——这是门 1 需要面对的诚实结论；'
    '若平滑窗口下尾巴出现与窗口无关的指数衰减 e^{-t/Δ}，则是真实修正。'
)

print(json.dumps(results, ensure_ascii=True, indent=2))
with open(Path(__file__).with_name('exp_cutoff_robust_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
