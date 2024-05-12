using Robin2011_RepPackage
using CSV
using DataFrames
using Distributions
using Random
using Optimization
using OptimizationNLopt
using DelimitedFiles

# Initialisation
global N = 150
global M = 500

## External Parameters
global r = 0.05/4;
global α = 0.5;
global τ = 0.5;
global k = 0.12;

dta = CSV.read("/Users/lewinnolden/Computational Economics/term_project/Robin 2011 Rep Files/USquarterly.csv", DataFrame, header = false)

filter!(row -> !isnan(row.Column1) && !isnan(row.Column4), dta);

prod = collect(dta.Column1);
prod = prod/mean(prod);
T = length(prod);
wage = collect(dta.Column2);
unr = collect(dta.Column3);
vac = collect(dta.Column4);
tightness = collect(dta.Column5);
f_q = collect(dta.Column6);
s_q = collect(dta.Column7);

X = [ prod unr s_q f_q ]

Me = mean(X, dims = 1)
Sd = std(log.(X), dims = 1)
kurt = [kurtosis(log.(X[:, j])) + 3 for i in 1:1, j in 1:4]
Co = cor(X)
elas_x = Co[1, :]'.*(Sd./Sd[1])
xreg = [ones(T-1) log.(prod[1:T-1])]
rhox = (xreg'xreg)\(xreg'*log.(prod[2:T]))

b0 = [Me[1], Sd[1], rhox[2], Me[2], Sd[2], kurt[2], Me[4]]

logit = x -> 1/(1 + exp(x))

b = parse.(Float64, readlines("x0.txt"))

ν = exp(b[1])
μ = exp(b[2])
δ = logit(b[3])
λ0 = logit(b[4])
ρ = tanh(b[5])
σ = exp(b[6])
z0 = exp(b[7])
C = exp(b[8])
B = 1

params_in = Dict(
    :ν => ν,
    :μ => μ,
    :δ => δ,
    :λ0 => λ0,
    :ρ => ρ,
    :σ => σ,
    :z0 => z0,
    :C => C
)

T = 5000;
burn = 1000;
Random.seed!(42)
draw = rand(burn+T, 1)

## Estimation

loginv = x -> log(1/x - 1)

x0 = [log(ν), log(μ), loginv(δ), loginv(λ0), atanh(ρ), log(σ), log(z0), log(C)]

f = OptimizationFunction((b, _) -> estCrit(b; draw = draw, burn = burn, b0 = b0, T = T, τ = τ, α = α, r = r, N = N, M = M))

prob = OptimizationProblem(f, x0)

sol = solve(prob, NLopt.LN_NELDERMEAD(); maxiters = 100000,reltol = 1e-4 )

params_opt = sol.u

params_opt2 = Dict(
    :ν  => exp(params_opt[1]),
    :μ => exp(params_opt[2]),
    :δ => logit(params_opt[3]),
    :λ0 => logit(params_opt[4]),
    :ρ => tanh(params_opt[5]),
    :σ => exp(params_opt[6]),
    :z0 => exp(params_opt[7]),
    :C => exp(params_opt[8])
)

writedlm("x1.txt", params_opt)