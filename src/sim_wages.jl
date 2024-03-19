using CSV
using DataFrames
using Distributions
using Robin2011_RepPackage

## Parameter Values - Change for when we estimate
δ::Number = 0.041563759920623
λ0::Number = 0.994544861919718
λ1::Number = 0.119345383430366
ρ = 0.913702234286476
ν::Number = 2.019365636076711
μ::Number = 5.786082109731152



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
year = collect(dta.Column8);
qtr = collect(dta.Column9);

time = year + qtr./4;

### Define Grids
include("grids.jl")

##Production and Surplus
p = matchprod(x, y)
z = homeprod(x, y)
Sx = SurplusVFI(p, z, Π)

## Wages
wd = WageVFI(Sx, Π, z)

##Steady State
ux = (δ/(δ + λ0)) .* (Sx .> 0) + (Sx .<= 0)

## Initial Conditions
burn = 0;
T1 = burn+T;
Ft = ones(T1)
uxt = ones(T1+1, M)
yt = ones(T1)
statet=zeros(T1)
Wqua=zeros(T1,9);
wqua=zeros(T1,9);
# wminqua=zeros(T1,9);
# wmaxqua=zeros(T1,9);
wagemean=zeros(T1);
wagevar=zeros(T1);
meanwvar=zeros(T1);
varwmean=zeros(T1);

i = min(sum(y.-prod[1] .<= 0), N);
Ft[1] = F[i];
yt[1] = y[i];
statet[1] = y[i];
uxt[1, :] = ux[i,:]

## Productivity and Unemployment Dynamics
for t in 1:(burn+T)

    uxt_1  = 1 .- (Sx .> 0) .* repeat((1 - δ) * (1 .- uxt[t, :]) + λ0 * uxt[t, :], 1, N)'
    ext = repeat(l', N, 1) .* (1 .- uxt_1)
    et = sum(ext, dims = 2)
    xt = vec(sum(p .* ext, dims = 2)./et)
    
    ut = uxt_1 * l
    
    i = argmin(abs.(xt .- prod[1]))
    statet[t] = i
    yt[t] = y[i]

    uxt[t+1, :] = [1 - (Sx[i, m] > 0) * ((1 - δ) * (1 - uxt[t, m]) + λ0 * uxt[t, m]) for m in 1:M]
end

#wdp = wage_dens_path(S, uxt, wd, l, U, statet, T)


