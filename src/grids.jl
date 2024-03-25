using Distributions
using Copulas

trim = 0.002
M = 500
N = 100
σ = 0.0257

x = collect(LinRange(trim, 1 - trim, M))
F = collect(LinRange(trim, 1- trim, N))
y = exp.(σ * quantile(Normal(), F))
l = pdf(Beta(ν, μ), x)
l = l ./ sum(l)

Σ = [1 ρ
     ρ 1]

cop = GaussianCopula(Σ)

P = [pdf(cop, [F[i], F[j]]) for i in 1:N, j in 1:N ]

Π = [P[i, j]/sum(P[i, :]) for i in 1:N, j in 1:N]
