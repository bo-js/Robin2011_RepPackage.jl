"""
`params_estimated(; path = "x1.txt")`

This function outputs all necessary parameters for the model, from a set of estimated parameters.

As it's input it takes a file path, which must point to a `.txt` file containing the parameter values correctly ordered, line separated, and transformed as in `estmain.jl`. 

The function then performs the appropriate inverse transformations and returns the parameters as a `NamedTuple`. 

The default file path points to those estimated by `estmain.jl`, stored in `x1.txt`. Alternatively to use those from Robin(2011), replace this with `x0.txt`. 

The parameters must be stored in the text file in the following order with the following transformations.

- `log(ν)`
- `log(μ)`
- `loginv(δ)`
- `loginv(λ0)`
- `atanh(ρ)`
- `log(σ)`
- `log(z0)`
- `log(C)`

Where `logit = x -> 1/(1+exp(x))` and `loginv = x -> log(1/x - 1)`.
"""
function params_estimated(; path = "x1.txt")
    b = parse.(Float64, readlines(path))
    
    logit = x -> 1/(1 + exp(x))

    ν = exp(b[1])
    μ = exp(b[2])
    δ = logit(b[3])
    λ0 = logit(b[4])
    ρ = tanh(b[5])
    σ = exp(b[6])
    z0 = exp(b[7])
    C = exp(b[8])

    r = 0.05/4
    τ = 0.5
    β = (1-δ)/(1+r)
    α = 0.64
    λ1 = 0.12 * λ0
    B = 1

    return (
        r = r,
        δ = δ,
        λ0 = λ0,
        λ1 = λ1,
        ρ = ρ,
        σ = σ,
        ν = ν,
        μ = μ,
        τ = τ,
        β = β,
        α = α,
        z0 = z0,
        C = C,
        B = B
    )
end
