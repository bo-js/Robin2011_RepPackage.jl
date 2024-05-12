
"""
`function unemp_path(S::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, δ::Number = 0.041563759920623)`

This function calculates the unemployment path for a given evolution of the aggregate state.

The Function takes as arguments:

- S, the surplus matrix, calculated using [`SurplusVFI`](@ref)
- statet, a given evolution of the aggregate state index through time
- T, the number of periods

As well as parameters:
- λ0, the rate at which an unemployed worker meets a firm 
- δ, the exogeneous job offer separatioin rate
The default parameter values are those used in Robin(2011).

The function returns a matrix, where each row gives a time period, and the columns give the unemployment rate for each worker type.
"""
function unemp_path(S::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    uxt = zeros(T, M)

    ## First Period Starts at Steady State
    uxt[1, :] = [(δ/(δ + λ0))*(S[statet[t], m] > 0) + 1 - (S[statet[t], m] > 0) for m in 1:M]

    ## Define Further Periods Recursively
    uxt[2:T, :] = [1 - ((1 - δ)*(1 - uxt[t-1, m]) + λ0 * uxt[t-1, m])*(S[statet[t-1], m] > 0) for t in 2:T, m in 1:M]

    return uxt

end



"""
`function wage_dens_path(S::Matrix, uxt::Matrix, wd::Dict, l::Vector, U::Matrix, statet::Vector, T::Integer; λ0::Number = 0.9945, λ1::Number = 0.1193, δ::Number = 0.0416)`

Wages are assigned through Bertrand Competition, resulting either in workers being paid their 
reservation wage (the monopsony wage), when being hired from unemployment, or capturing the full surplus when being poached. 
However, since workers do not receive a poaching offer every period, you may observe multiple wages at every single point in time.
Furthermore, changes in the state may make it so that either the worker or the firm have a credible threat to leave/ fire leading to renegotiation of the wage.

This function calculates the wage density path for a given evolution of the aggregate state, according to the laws of motion derived in Robin (2011),
giving the density of workers earning each of the possible wages at each point in time.

The Function takes as arguments:

- `S`, the surplus matrix, calculated using [`SurplusVFI`](@ref)
- `uxt`, the unemployment path, calculated using [`unemp_path'](@ref)
- `wd`, a Dict including the minimum and the maximum wage in each period as well as their current value, calculated using [`WageVFI`](@ref)
- `l`, the number of workers of each type, calculated using [`grids`](@ref)
- `U`, the value of unemployment 
- `statet`, a given evolution of the aggregate state index through time
- `T`, the number of periods

As well as parameters:

- `λ0`, the rate at which an unemployed worker meets a firm 

The default parameter values are those used in Robin(2011).

The function returns a four dimensional array, where each entry is the measure of workers earning that wage. 
The first index denotes the time period, the second the state that wage is assigned, the third the worker type that wage is assigned to, and the fourth wether the wage is the monopsony wage `[:,:,:1]`,
or the full surplus wage `[:, :, :, 2]`.
"""
function wage_dens_path(S::Matrix, uxt::Matrix, wd::Dict, l::Vector, U::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, λ1::Number = 0.119345383430366, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    N = length(S[:, 1])

    gt = zeros(T, N, M, 2)

    Wmin = wd[:Wmin]
    Wmax = wd[:Wmax]

    # Define First Period from Steady State
    gt[1, statet[1], :, 1] = [(λ0 * (S[statet[1], m] > 0) * uxt[1, m] * l[m])/(1 - (1 - δ) * (S[statet[1], m] > 0) * (1 - λ1)) for m in 1:M]
    
    gt[1, statet[1], :, 2] = [((1 - δ) * λ1 * (1 - uxt[1, m]) * (S[statet[1], m] > 0) * l[m])/(1 - (1 - δ) * (S[statet[1], m] > 0) * (1 - λ1)) for m in 1:M]

    # Define subsequent Periods Recurively
    for t in 2:T
        for i in 1:N
            if statet[t] == i # For wages being assigned in current period's state
                
                gt[t, i, :, 1] = [(S[statet[t], m] > 0) * (λ0 * uxt[t-1, m] * l[m] + (1 - δ) * (1 - λ1) * (
                    gt[t-1, i, m, 1] + sum(
                        
                        (Wmin[statet[t], k, m] - U[statet[t], m] < 0) * gt[t-1, k, m, 1] + 
                        
                        (Wmax[statet[t], k, m] - U[statet[t], m] < 0) * gt[t-1, k, m, 2] for k in 1:N
                    )
                )) for m in 1:M]

                gt[t, i, :, 2] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (
                        λ1 * (1 - uxt[t - 1, m]) * l[m] + (1 - λ1) * (
                            gt[t - 1, i, m, 2] + sum(
                                
                                (Wmin[statet[t], k, m] - U[statet[t], m] > S[statet[t], m]) * gt[t-1, k, m, 1] + 
                        
                                (Wmax[statet[t], k, m] - U[statet[t], m] > S[statet[t], m]) * gt[t-1, k, m, 2] for k in 1:N
                            )
                        )
                    )
                for m in 1:M]

            else # For wages still viable in this state, but not being actively assigned.

                gt[t, i, :, 1] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (1 - λ1) * (0 ≤ Wmin[statet[t], i, m] - U[statet[t], m] ≤ S[statet[t], m])
                    for m in 1:M
                ]

                gt[t, i, :, 2] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (1 - λ1) * (0 ≤ Wmax[statet[t], i, m] - U[statet[t], m] ≤ S[statet[t], m])
                    for m in 1:M
                ]

            end
        end
    end

    return gt

end