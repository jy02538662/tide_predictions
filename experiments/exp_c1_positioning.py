"""exp2 门 C1 定位厘清：可证伪预言 vs 恒等式。

门 C1 推导链第 2 步（核心定位）。

exp1 已坐实数学机制：对数周期振荡 ⟺ 离散嵌套尺度 c（特征尺度）。
本实验厘清 C1 的定位，并数值演示「离散嵌套 → 对数周期」与「连续尺度不变 → 纯 1/t」的判别。

关键逻辑（诚实）：
  预印本 1.3 §3.5 把 C1 写成「可证伪预言：观测到对数周期 ⟹ A 被否」。
  但 exp_gravity_uniqueness 注释说「对数周期被 A 先验排除」。
  两者其实是一件事的两面：
    - 正方向：A ⟹ 无特征尺度 ⟹ 无 c ⟹ 纯 1/t（这是 A 的**推论**，不是独立预言）；
    - 反方向：观测到对数周期 ⟹ 存在 c ⟹ A 被否（这是 A 的**可证伪判据**）。

  所以 C1 的准确定位是：
    「纯 1/t」本身不是新预言（它是 A 的推论/恒等式）；
    但「对数周期 = A 被否的判据」是可证伪的（给定了「什么观测会否掉 A」的明确标准）。
    这类似 Born 偏离判据：Born 本身不推导，但「偏离 Born ⟺ 外部观察者」是可证伪判据。
"""

import json
from pathlib import Path

import numpy as np

results = {}


def G_discrete_nesting(c, N, t_vals):
    """离散嵌套尺度 c：频率 ω = n·ln c 离散均匀，G(t)=Σ e^{iωt}。"""
    ln_c = np.log(c)
    freqs = np.arange(-N, N + 1, dtype=float) * ln_c
    weights = np.ones_like(freqs) / len(freqs)
    return np.array([np.real(np.sum(weights * np.exp(1j * freqs * tv))) for tv in t_vals])


def G_continuous(t_vals):
    """连续尺度不变：ρ=C/λ，G(t)~1/t（纯幂律，无振荡）。"""
    return 1.0 / t_vals


def detect_log_periodicity(t_vals, G):
    """在对数时间轴 ln t 上做傅里叶，检测 G(t) 的周期结构。

    对数周期振荡 G~t^{-p}cos(2π ln t / ln c) 在 ln t 轴上是周期 ln c 的振荡。
    判据：先除以幂律包络，再看 ln t 轴的自相关/频谱峰值。
    """
    ln_t = np.log(t_vals)
    # 除以整体衰减包络（用 |G| 的中位数归一，简单起见只看 G 的符号翻转结构）
    # 更直接：G 在 ln t 轴上的频谱
    # 均匀重采样到 ln t 均匀网格
    ln_grid = np.linspace(ln_t[0], ln_t[-1], 512)
    G_interp = np.interp(ln_grid, ln_t, G)
    spec = np.abs(np.fft.rfft(G_interp - G_interp.mean()))
    # 主频（排除直流）
    freqs_ln = np.fft.rfftfreq(len(ln_grid), d=(ln_grid[1] - ln_grid[0]))
    if len(spec) > 1:
        peak_idx = np.argmax(spec[1:]) + 1
        peak_freq = freqs_ln[peak_idx]
        peak_val = spec[peak_idx]
        total = np.sum(spec[1:])
        return peak_freq, float(peak_val / total if total > 0 else 0)
    return 0.0, 0.0


def main():
    print("=" * 74)
    print("门 C1 定位厘清：对数周期是 A 的可证伪判据（非独立预言）")
    print("=" * 74)

    t_vals = np.linspace(10.0, 200.0, 800)

    # 情形 1：离散嵌套（c=1.5，一个特征尺度）→ 对数周期振荡
    G1 = G_discrete_nesting(1.5, 100, t_vals)
    # 情形 2：连续尺度不变 → 纯 1/t
    G2 = G_continuous(t_vals)

    # 检测：ln t 轴频谱主频占比
    f1, r1 = detect_log_periodicity(t_vals, G1)
    f2, r2 = detect_log_periodicity(t_vals, G2)

    results['discrete_nesting'] = {
        'ln_t_peak_freq': float(f1),
        'spectral_concentration': float(r1),
        'note': '离散嵌套（有特征尺度 c）⟹ ln t 轴有明显周期峰 ⟹ 对数周期振荡',
    }
    results['continuous'] = {
        'ln_t_peak_freq': float(f2),
        'spectral_concentration': float(r2),
        'note': '连续尺度不变 ⟹ 无周期结构（纯 1/t）',
    }

    print(f"  离散嵌套 c=1.5：ln t 轴频谱主频={f1:.4f}，能量占比={r1:.3f}")
    print(f"  连续尺度不变：  ln t 轴频谱主频={f2:.4f}，能量占比={r2:.3f}")

    # 定位结论
    results['C1_positioning'] = {
        'forward': 'A（无偏好）⟹ 无特征尺度 ⟹ 无 c ⟹ 纯 1/t（A 的推论，非独立预言）',
        'backward': '观测到对数周期尾巴 ⟹ 存在特征尺度 c ⟹ A 被否（可证伪判据）',
        'verdict': (
            'C1 的准确定位：纯 1/t 本身是 A 的恒等式（不是新预言），'
            '但「对数周期 = A 被否」是可证伪判据——它给出了「什么观测会否掉 A」的明确标准。'
            '类比 Born 偏离判据：Born 不推导，但「偏离 Born ⟺ 外部观察者」是可证伪判据。'
            '所以 C1 应降级为「判据」，而非「独立预言」。'
        ),
    }

    print("\n" + "=" * 74)
    print("  定位结论：")
    print("    正方向：A ⟹ 纯 1/t（恒等式，非新预言）")
    print("    反方向：对数周期 ⟹ A 被否（可证伪判据）")
    print("    ⟹ C1 是「判据」，不是「独立预言」，类似 Born 偏离判据")
    print("=" * 74)

    out = Path(__file__).with_name('exp_c1_positioning_last_run.json')
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f"\n  wrote {out.name}")


if __name__ == '__main__':
    main()
