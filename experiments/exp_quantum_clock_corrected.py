"""门 4 量子时钟修正：诚实厘清——「时间分辨率 ∝ 模流谱隙」的猜测不成立。

门 4 推导链（诚实版，纠正猜测）。

v9 的猜测：「时间 = 模流 σ_t ⟹ 时间分辨率依赖模流谱」。

本实验发现这个猜测混淆了两个不同的量：
  - 模流谱的「谱隙」Δθ = 最小 |θ_ij|（最小频率间隔）→ 决定**最慢振荡（准周期 T=2π/Δθ）**；
  - 模流谱的「谱宽」= ln(λ_max/λ_min)（频率分布范围）→ 决定**最快振荡**。

时间分辨率（能测多小的时间差）由**能量-时间不确定关系 δt ~ 1/δE** 决定，
δE 是**能量展宽（谱宽）**，不是谱隙。所以「时间分辨率 ∝ 谱隙」的猜测
物理上站不住——谱隙对应「准周期（最慢振荡）」，不是「时间分辨率」。

结论：门 4 的「模流谱 → 时间分辨率」猜测**方向有误**（谱隙 ≠ 谱宽）。
「时间=模流」是框架独有结构（这个没错），但「时间分辨率依赖谱隙」
这条具体猜测链走不通。
"""

import json
from pathlib import Path

import numpy as np

results = {}

# 数值坐实：谱隙 vs 谱宽的区分
N = 64
rng = np.random.default_rng(0)
lam = np.sort(rng.uniform(0.5, 2.0, N))
theta = np.log(lam[:, None]) - np.log(lam[None, :])
theta_pos = np.abs(theta[theta > 0])
spectral_width = float(np.max(theta))       # ln(λ_max/λ_min)
spectral_gap = float(np.min(theta_pos))     # 最小频率间隔

results['gap_vs_width'] = {
    'spectral_width_ln_lmax_lmin': spectral_width,
    'spectral_gap_min_theta': spectral_gap,
    'note': (
        f'谱宽 = ln(λ_max/λ_min) = {spectral_width:.4f}（决定最快振荡，时间分辨率由它定）；'
        f'谱隙 = 最小|θ| = {spectral_gap:.2e}（决定最慢振荡 = 准周期，不是时间分辨率）。'
    ),
}

results['corrected_physics'] = {
    'time_resolution': 'δt ~ 1/δE，δE = 能量展宽（谱宽），来自能量-时间不确定关系',
    'spectral_gap_meaning': '谱隙 Δθ → 准周期 T=2π/Δθ（最慢振荡），不是时间分辨率',
    'confusion': 'v9 猜测把「谱隙」当「谱宽」，故「时间分辨率 ∝ 谱隙」不成立',
}

results['verdict'] = {
    'conclusion': (
        '门 4 = 猜测方向有误（负结果）。「时间=模流」是框架独有结构（保留），'
        '但「时间分辨率依赖模流谱隙」这条具体猜测链走不通（谱隙 ≠ 谱宽）。'
    ),
    'what_survives': (
        '框架独有的是「时间 = 模流（从观察者态涌现）」，不是「时间分辨率依赖谱隙」。'
        '「时间分辨率」若要成为预言，需另找「模流谱宽 → 时间分辨率」的正确链条，'
        '但那是能量-时间不确定关系的重述，不独有。'
    ),
    'strength': '负结果（比门 1 更弱——门 1 是「结构上无修正」，门 4 是「猜测本身混淆了概念」）',
}

print(json.dumps(results, ensure_ascii=False, indent=2))
with open(Path(__file__).with_name('exp_quantum_clock_corrected_last_run.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
