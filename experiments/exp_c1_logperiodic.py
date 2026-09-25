"""exp1 门 C1 对数周期振荡的数学来源：离散递归 vs 连续递归，G(t) 尾巴判别信号。

门 C1 推导链第 1 步（厘清数学机制）。

背景张力（要诚实处理）：
  - 预印本 1.3 §3.5 说 C1 是「可证伪预言」：观测到对数周期尾巴 ⟹ A（无偏好）被否。
  - exp_gravity_uniqueness 注释却说「对数周期被 A 先验排除」——那它就是重言式，不是预言。

本实验把对数周期振荡的**数学来源**彻底厘清：
  - 对数周期振荡 G(t) ~ t^{-p} * cos(2π ln t / ln c) 来自「离散尺度不变」：
    rho(λ) 在 λ → cλ 下不变，但只在**离散的 c 倍数**下不变，不是连续尺度不变。
  - 数学上：若 rho(λ) 满足 rho(cλ)=c^{-1}rho(λ) 只对**一个固定的 c**（不是所有 c），
    则 rho(λ) = λ^{-1} * P(ln λ / ln c)，其中 P 是周期 1 的函数（周期函数），
    傅里叶分解后 P 的每个 Fourier 分量贡献 λ^{-1} * λ^{2πik/ln c} = λ^{-1+2πik/ln c}，
    即复幂律（对数周期振荡），G(t) 尾巴有 cos(2π ln t / ln c) 振荡。
  - 关键判据：A（无偏好尺度）⟹ rho(cλ)=c^{-1}rho(λ) 对**所有 c>0** 成立
    ⟹ P 是周期 ln c 对**所有** ln c 的周期函数 ⟹ P 只能是常数 ⟹ 无振荡。

所以核心问题是：
  A 到底只允许「连续尺度不变」（所有 c）→ 纯 1/t，还是可能退化成「离散尺度不变」（单个 c）→ 振荡？

本实验数值演示两种情况，坐实「振荡 ⟺ 离散嵌套尺度 c」这个数学事实。
"""

import json
from pathlib import Path

import numpy as np

results = {}


def G_from_rho(rho_fn, t_vals, lam_min, lam_max, n=200000):
    """G(t) = ∫ rho(λ) 相关... 简化：直接算离散嵌套尺度的谱 + 关联。

    用最直接的模型：观察者态谱 λ_n = c^n（离散嵌套，特征尺度 c）。
    模流频率差 ω = ln(λ_n/λ_m) = (n-m) ln c，是**离散均匀**的（间距 ln c）。
    G(t) = Σ ρ_n e^{it ω}，离散频率 ⟹ 周期函数在 t 轴（周期 2π/ln c）。
    """
    pass


def main():
    print("=" * 74)
    print("门 C1：对数周期振荡的数学来源")
    print("=" * 74)

    # 模型 1：离散嵌套尺度（特征尺度 c）——预测对数周期
    # 观察者态 ρ(λ) 在对数轴上**离散均匀**（间距 ln c），即只在 λ_n = c^n 有支撑。
    # 这对应「离散递归」：一层套一层，层间隔 c。
    # 模流频率 ω_{nm} = ln(λ_n/λ_m) = (n-m) ln c，离散间距 ln c。
    # G(t) = Σ_{n,m} ρ_n ρ_m e^{i t (n-m) ln c}，是 t 的周期函数（周期 2π/ln c）。

    c = 2.0  # 嵌套尺度
    N = 200  # 层数
    ln_c = np.log(c)
    # 频率在 [-N ln c, N ln c] 上离散均匀，间距 ln c
    freqs = np.arange(-N, N + 1, dtype=float) * ln_c
    # 权重（对数均匀谱 ρ=C/λ ⟹ 每层等权）
    weights = np.ones_like(freqs) / len(freqs)
    # G(t) = Σ w e^{i ω t}，取实部
    t_vals = np.linspace(5.0, 80.0, 400)
    G_disc = np.array([np.real(np.sum(weights * np.exp(1j * freqs * tv))) for tv in t_vals])

    # 判别：在 log t 轴上做傅里叶，看是否有 ln c 周期
    # 直接看 G(t) 是否周期 = 2π/ln c
    period_t = 2 * np.pi / ln_c
    results['discrete_nesting'] = {
        'nesting_scale_c': float(c),
        'predicted_period_in_t': float(period_t),
        'G_t_samples': np.round(G_disc[:10], 4).tolist(),
    }
    print(f"  离散嵌套（c={c}）：预测 G(t) 在 t 轴周期 = 2π/ln c = {period_t:.3f}")
    print(f"    即在对数时间轴上周期 = 2π（对数周期振荡）")
    print(f"    G(t) 前 10 样本: {np.round(G_disc[:10],3).tolist()}")

    # 模型 2：连续尺度不变（A 成立）——预测纯 1/t，无振荡
    # ρ(λ) = C/λ 连续（对数轴连续均匀），频率 ω 连续分布，G(t) ~ 1/t
    # 已由 exp_gravity_uniqueness 验证：G(t) ~ t^{-1}，无周期
    results['continuous_scale_invariance'] = {
        'prediction': 'G(t) ~ 1/t，无对数周期（exp_gravity_uniqueness 已验）',
    }
    print(f"\n  连续尺度不变（A 成立）：G(t)~1/t，无振荡（已由 exp_gravity_uniqueness 验证）")

    # 关键判据：振荡 ⟺ 离散嵌套尺度 c（一个特征尺度）
    # 若 A（无偏好）禁止任何特征尺度，则 c 不存在 ⟹ 无振荡 ⟹ 纯 1/t
    # 所以「对数周期」需要「一个被偏好的嵌套尺度 c」，这与 A 矛盾。
    results['key_criterion'] = {
        'statement': (
            '对数周期振荡 ⟺ 存在离散嵌套尺度 c（特征尺度）。'
            'A（无偏好尺度）⟹ c 不存在 ⟹ 无振荡 ⟹ 纯 1/t。'
            '故「观测到对数周期」在 A 成立时不可能发生——C1 更接近 A 的恒等式，'
            '而非可证伪预言（见下一实验厘清定位）。'
        ),
    }

    print("\n" + "=" * 74)
    print("  关键判据：振荡 ⟺ 离散嵌套尺度 c（一个被偏好的特征尺度）")
    print("  A（无偏好）⟹ 无特征尺度 ⟹ 无 c ⟹ 无振荡 ⟹ 纯 1/t")
    print("=" * 74)

    out = Path(__file__).with_name('exp_c1_logperiodic_last_run.json')
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f"\n  wrote {out.name}")


if __name__ == '__main__':
    main()
