function WageVFI(S::Matrix, Π::Matrix, z::Matrix; λ1::Number = 0.119345383430366, β::Number = 0.946603693905558)
    
    N = length(Π[1, :])
    M = length(S[1, :])

    Wmin = zeros(N, N, M)
    Wmax = zeros(N, N, M)
    wmin = zeros(N, M)
    wmax = zeros(N, M)

    for m in 1:length(M)
        for i in 1:length(N)

            # Minimum Wages - W is a vector giving the worker surplus of being in each state with the minimum wage set for state i
            A = (Π - repeat(Π[i, :]', N, 1)) * repeat((S[:, m] .> 0)', N, 1)
            W = (I(length(Π[1, :])) - β * A * (1 - λ1)) \ (z[i, m] .- z[:, m] + β * A * λ1 * S[:, m])
            W0 = min.(max.(W, 0), S[:, m])
            e = norm(W - W0, 2)

            while e > 0.00001
                W1 = W
                W = z[i, m] .- z[:, m] + β * A * (λ1 * S[:, m] + (1 - λ1) * W0)
                W0 = min.(max.(W, 0), S[:, m])
                e = norm(W - W1, 2)
            end

            Wmin[:, i, m] = W
            wmin[i, m] = z[i, m] - β * Π[i, :]'*(S[:, m] .> 0)


            ## Maximum Wages - Likewise as above

            W = (I(length(Π[1, :])) - β * A * (1 - λ1)) \ (S[i, m] + z[i, m] .- z[:, m] + β * A * λ1 * S[:, m])
            W0 = min.(max.(W, 0), S[:, m])
            e = norm(W - W0, 2)

            while e > 0.00001
                W1 = W
                W = S[i, m] + z[i, m] .- z[:, m] + β * A * (λ1 * S[:, m] + (1 - λ1) * W0)
                W0 = min.(max.(W, 0), S[:, m])
                e = norm(W - W1, 2)
            end

            Wmax[:, i, m] = W
            wmax[i, m] = S[i, m] + z[i, m] - β * Π[i, :]'*(S[:, m] .> 0)


        end
    end

    return Dict(:Wmin => Wmin, :wmin => wmin, :Wmax => Wmax, :wmax => wmax)
end