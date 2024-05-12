"""
`WageVFI(S::Matrix, Π::Matrix, z::Matrix; λ1 = 0.119, β::Number = 0.9466)`

This function calculates equilibrium wages, as well as their associated value functions.

Wages are assigned through Bertrand Competition, resulting either in workers being paid their reservaiton wage (the monopsony wage), when being hired from unemployment, or capturing the full surplus when being poached.

The Function returns a Dict with keys:
- `:wmin`, an `i x j` matrix giving the monopsony wage assigned to worker `j`, being hired out of unemployment in state `i`.
- `:wmax`, an `i x j` matrix giving the full surplus wage assigned to worker `j`, being poached in state `i`.
- `:Wmin`, an `i x j x k` array giving the value of worker `k`, being in state `i`, with the monopsony wage from state `j`.
- `:Wmax`, an `i x j x k` array giving the value of worker `k`, being in state `i`, with the poacher's wage from state `j`.

The Function takes as arguments:
- `S`, the surplus matrix, calculated using [`SurplusVFI`](@ref).
- `Π`, the markov transition matrix, calculated using [`grids`](@ref).
- `z`, the home production matrix, calculated using [`homeprod`](@ref).

As well as parameters:
- `λ1`, the probability of meeting a new employer while already employed.
- `β`, the discounted probability of exogenenous separation, i.e. `(1 - δ)/(1+r)`.

The default parameter values are those used in Robin(2011).
"""
function WageVFI(S::Matrix, Π::Matrix, z::Matrix; λ1 = 0.119345383430366, β::Number = 0.946603693905558)
    
    N = length(Π[1, :])
    M = length(S[1, :])

    Wmin = zeros(N, N, M)
    Wmax = zeros(N, N, M)
    wmin = zeros(N, M)
    wmax = zeros(N, M)
    kf = zeros(N)
    kf .= λ1
    for m in 1:M
        for i in 1:N

            # Minimum Wages - W is a vector giving the worker surplus of being in each state with the minimum wage set for state i
            A = (Π - repeat(Π[i, :]', N, 1)) * repeat((S[:, m] .> 0)', N, 1)
            W = (I(length(Π[1, :])) - β * A * repeat((1 .- kf)', N, 1)) \ (z[i, m] .- z[:, m] + β * A *( kf .* S[:, m]))
            W0 = min.(max.(W, 0), S[:, m])
            e = norm(W - W0, 2)

            while e > 0.001
                W1 = W
                W = z[i, m] .- z[:, m] + β * A * (kf .* S[:, m] + (1 .- kf) .* W0)
                W0 = min.(max.(W, 0), S[:, m])
                e = norm(W - W1, 2)
            end

            Wmin[:, i, m] = W
            wmin[i, m] = z[i, m] - first(β * Π[i, :]' .* (S[:, m] .> 0)' * (kf .* S[:, m] + (1 .- kf) .* W0))


            ## Maximum Wages - Likewise as above

            W = (I(length(Π[1, :])) - β * A * repeat((1 .- kf)', N, 1)) \ (S[i, m] + z[i, m] .- z[:, m] + β * A * (kf .* S[:, m]))
            W0 = min.(max.(W, 0), S[:, m])
            e = norm(W - W0, 2)

            while e > 0.001
                W1 = W
                W = S[i, m] + z[i, m] .- z[:, m] + β * A * (kf .* S[:, m] + (1 .- kf) .* W0)
                W0 = min.(max.(W, 0), S[:, m])
                e = norm(W - W1, 2)
            end

            Wmax[:, i, m] = W
            wmax[i, m] = S[i, m] + z[i, m] - first(β * Π[i, :]' .* (S[:, m] .> 0)' * (kf .* S[:, m] + (1 .- kf) .* W0))


        end
    end
 
    return Dict(:Wmin => Wmin, :wmin => wmin, :Wmax => Wmax, :wmax => wmax)
end