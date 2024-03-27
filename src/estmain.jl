using Robin2011_RepPackage
using CSV
using DataFrames
using Distributions

# Initialisation
global N = 150
global M = 500

## External Parameters
global r = 0.05/4;
global α = 0.5;
global τ = 0.5;
global k = 0.12;

dta = CSV.read("/Users/bojs/Desktop/Robin 2011 Rep Files/matlab/USquarterly.csv", DataFrame, header = false)

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
λ0 = logit(b[3])
ρ = tanh(b[5])
σ = exp(b[6])
z0 = exp(b[7])
C = exp(b[8])
B = 1

T = 5000;
burn = 1000;
draw = rand(burn+T, 1)

## Estimation

loginv = x -> log(1/x - 1)

x0 = [log(ν), log(μ), loginv(δ), loginv(λ0), atanh(ρ), log(σ), log(z0), log(C)]

estCrit(b; burn = burn)