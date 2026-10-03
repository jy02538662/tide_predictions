import sympy as sp
from sympy import symbols, Function, diff, simplify, sqrt, Matrix, Rational

x = symbols('x')
c = Function('c')(x)   # 声速场 c(x)，c² ∝ 密度 ρ
v = Function('v')(x)   # 流速场 v(x)

# --- 1+1D 声学度规（Painlevé-Gullstrand 型）：ds² = -(c²-v²)dt² - 2v dx dt + dx² ---
# 度规张量 g_μν，指标 (t, x)
g = Matrix([
    [-(c**2 - v**2), -v],
    [-v, 1],
])
ginv = g.inv()

# Christoffel Γ^ρ_{μν} = 1/2 g^{ρσ}(∂_μ g_{νσ} + ∂_ν g_{μσ} - ∂_σ g_{μν})
# 2D 指标 (t=0, x=1)
G = [[[0,0],[0,0]], [[0,0],[0,0]]]  # G[ρ][μ][ν]

def gd(mu, nu):
    # ∂_μ g_{νσ} 等：只有 x 导数非零
    d = sp.diff(g[nu, mu], x) if mu == 1 else 0
    return d

for rho in range(2):
    for mu in range(2):
        for nu in range(2):
            s = 0
            for sig in range(2):
                s += ginv[rho, sig] * (gd(mu, nu) if False else 0)
            # 直接公式
            term = 0
            for sig in range(2):
                t = sp.diff(g[nu, sig], x) * (1 if mu==1 else 0) \
                  + sp.diff(g[mu, sig], x) * (1 if nu==1 else 0) \
                  - sp.diff(g[mu, nu], x) * (1 if sig==1 else 0)
                term += ginv[rho, sig] * t
            G[rho][mu][nu] = sp.simplify(sp.Rational(1,2) * term)

# Ricci 标量 R = g^{μν} R_{μν}, R_{μν} = ∂_ρ Γ^ρ_{μν} - ∂_ν Γ^ρ_{μρ} + ΓΓ - ΓΓ
R_mu = [[0,0],[0,0]]
for mu in range(2):
    for nu in range(2):
        val = 0
        for rho in range(2):
            val += sp.diff(G[rho][mu][nu], x) * (1 if rho==1 else 0)
            val -= sp.diff(G[rho][mu][rho], x) * (1 if nu==1 else 0)
        # ΓΓ 项（2D，写全）
        for rho in range(2):
            for lam in range(2):
                val += G[rho][mu][nu]*G[lam][rho][lam]
                val -= G[rho][mu][lam]*G[lam][nu][rho]
        R_mu[mu][nu] = sp.simplify(val)

R = sp.simplify(ginv[0,0]*R_mu[0][0] + ginv[0,1]*R_mu[0][1] + ginv[1,0]*R_mu[1][0] + ginv[1,1]*R_mu[1][1])
print("1+1D 声学度规标量曲率 R =", sp.simplify(R))

# 静态极限 v=0：R = ?
R_static = sp.simplify(R.subs(v, 0))
print("静态 v=0 极限 R =", R_static)
print("(预期 R = -2c''/c，即声速二阶导 / 声速)")
