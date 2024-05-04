function estCrit(b; M = 500, N = 100, τ = 0.5, α = 0.5, r = 0.05/4, T = 5000, burn = 1000, draw = rand(burn+T, 1), b0 = [1, 0.0226, 0.9136, 0.0579, 0.2141, 2.5296, 0.7842])

    logit = x -> 1/(1 + exp(x))

    ν = exp(b[1])
    μ = exp(b[2])
    δ = logit(b[3])
    λ0 = logit(b[4])
    ρ = tanh(b[5])
    σ = exp(b[6])
    z0 = exp(b[7])
    C = exp(b[8])
    B = 1

    g = grids(; M = M, N = N, ν = ν, μ = μ, ρ = ρ, σ = σ)
    x = g[:x]
    y = g[:y]
    Π = g[:Π]
    l = g[:l]

    p = matchprod(x, y; B = 1, C = C)

    z = homeprod(x, y; B = 1, C = C, α = 0.5, z0 = z0)

    Ux = (I(N) - Π./(1 + r))\z

    Sx = SurplusVFI(p, z, Π; β = (1 - δ)/(1+r))

    S = max.(Sx, 0)*l # Aggregate Surplus
    L = (Sx .> 0)*l

    ## Steady State
    ux = (δ/(δ + λ0)).*(Sx .> 0) + (Sx .≤ 0)
    u = 1 .- δ .* L ./ (δ + λ0)

    ## Dynamics

    i = Integer(round(N/2))
    yt = y[i] .* ones(T+burn)
    uxt = repeat(ux[i, :]', T+burn+1, 1)
    Sxt = repeat(Sx[i, :]', T+burn+1, 1)

    for t = 1:T+burn
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), N)
        yt[t] = y[i]
        Sxt[t, :] = Sx[i, :]
        uxt[t+1, :] = 1 .- (Sxt[t, :] .> 0) .* ((1 - δ) .* (1 .- uxt[t, :]) + λ0 .* uxt[t, :])
    end

    sel = (burn+1):(T+burn);
    yt = yt[sel]
    uxt = uxt[sel, :]
    Sxt = Sxt[sel, :]
    ut = uxt * l

    # Exit Rate from Unemployment
    ft = ((((Sxt .> 0).*uxt) * l) .* λ0) ./ut
    # Quit Rate
    qt = (((Sxt .> 0).*(1 .- uxt))*l) .* (0.12 * λ0 * τ) ./(1 .- ut)
    # Layoff Rate
    sxt = (1 .- (Sxt .> 0).*(1 - δ)).*(1 .- uxt)
    st = (sxt * l)./(1 .- ut)
    # Employment
    ext = repeat(l', T, 1).*(1 .- uxt)
    et = sum(ext, dims = 2)
    
    # Match Productivity
    pt = matchprod(x, yt; B = B, C = C)
    xt = sum(pt.*ext, dims = 2)./et
    
    # Leisure Cost
    zxt = homeprod(x, yt; B = B, C = C)
    zt = sum(Sxt.*ext, dims = 2)./et

    # Output

    X = [xt ut st ft]

    Me = mean(X, dims = 1)
    Sd = std(log.(X), dims = 1)
    kurt = [kurtosis(log.(X[:, j])) + 3 for i in 1:1, j in 1:4]
    Co = cor(log.(X))
    elas_x = Co[1, :]'.*(Sd./Sd[1])
    xreg = [ones(T-1) log.(xt[1:T-1])]
    rhox = (xreg'xreg)\(xreg'*log.(xt[2:T]))

    b1 = [Me[1], Sd[1], rhox[2], Me[2], Sd[2], kurt[2], Me[4]]

    db = (b1-b0)./b0
    crit = db'*db

    return crit
end