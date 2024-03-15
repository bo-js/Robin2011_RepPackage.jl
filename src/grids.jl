using Distributions
using BivariateCopulas

trim = 0.002
M = 500
N = 100
σ = 0.0257

x = collect(LinRange(trim, 1 - trim, M))
F = collect(LinRange(trim, 1- trim, N))
y = exp.(σ * quantile(Normal(), F))

ρ = 0.913702234286476

cop = Gaussian(ρ)

P = [density(cop, F[i], F[j]) for i in 1:N, j in 1:N ]

Π = [P[i, j]/sum(P[i, :]) for i in 1:N, j in 1:N]

μ = 5.786082109731152
ν = 2.019365636076711

l = pdf(Beta(μ,ν), x)
l = l./sum(l)





