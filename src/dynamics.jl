function unemp_path(S::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    ut = zeros(T, M)

    ## First Period Starts at Steady State
    ut[1, :] = [(δ/(δ + λ0))*(S[statet[t], m] > 0) + 1 - (S[statet[t], m] > 0) for m in 1:M]

    ## Define Further Periods Recursively
    ut[2:T, :] = [1 - ((1 - δ)*(1 - ut[t-1, m]) + λ0 * ut[t-1, m])*(S[statet[t-1], m] > 0) for t in 2:T, m in 1:M]

end