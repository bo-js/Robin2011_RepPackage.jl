
using Distributions
using BivariateCopulas

trim = 0.002
M = 500
N = 100

x = collect(LinRange(trim, 1-trim, M))

F = LinRange(trim, 1-trim, N)

# Manual meshgrid creation
U1, U2 = repeat(reshape(F, 1, :), length(F), 1), repeat(reshape(F, :, 1), 1, length(F))

sig = 0.0257

y = exp.(sig * quantile(Normal(), F))
ρ = 0.9137

cop = Gaussian(ρ)

P = [density(cop, F[i], F[j]) for i in 1:N, j in 1:N]

# Calculate the sum of each row in P
rowSums = sum(P, dims=2)

# Normalize P such that each row sums to 1
Π = P ./ rowSums

p = buildingblocks.matchprod(x,y)
z = buildingblocks.homeprod(x,y)
buildingblocks.SurplusVFI(p,z,Π) 