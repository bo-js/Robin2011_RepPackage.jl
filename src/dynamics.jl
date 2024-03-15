function unemp_path(S::Matrix, yt::Vector, T::Integer; λ0::Number = 0.994544861919718, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    ut = zeros(T, M)

    ut[1, :] = [(δ/(δ + λ0))*(S[yt[t], m] > 0) + 1 - (S[yt[t], m] > 0) for m in 1:M]

    ut[2:T, :] = [1 - ((1 - δ)*(1 - ut[t-1, m]) + λ0 * ut[t-1, m])*(S[yt[t-1], m] > 0) for t in 2:T, m in 1:M]

end