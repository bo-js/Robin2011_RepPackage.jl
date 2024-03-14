module Robin2011_RepPackage

using LinearAlgebra
using Distributions

function matchprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.7258)
    return y * (B * x .+ C)'
end

function homeprod(x::Vector, y::Vector;  B::Number = 1, C::Number = 0.7258, α::Number = 0.75, zθ::Number = 0.7668)
    return zθ .+ α * (matchprod(x,y) .- zθ)
end


function SurplusVFI(p::Matrix, z::Matrix, Π::Matrix; β::Number = 0.9466)
    S =(I(length(Π[1,:])) - (β * Π))\(p - z) 
    e = norm(S - max.(S,0),2)

    while e > 0.00001
        S1 = S
        S = p - z + β * Π * max.(S,0)
        e = norm(S - S1, 2)
    end

    return S 

end

function SteadyStateUnempl(Sx::Matrix; δ::Number = 0.041563759920623, ϕ::Number = 0.994544861919718)
    # use SurplusVFI to get Sx
    ux = (δ / (δ + ϕ)) * (Sx .> 0) + (Sx .<= 0)
end

Sx = SurplusVFI(p,z,Π)

ux = SteadyStateUnempl(Sx)

export matchprod
export homeprod
export SurplusVFI
export SteadyStateUnempl



# turnover dynamics

    # parameters for grid

    burn = 0
    T = 239
    T1 = burn + T
    uxt = ones(T1 + 1, M)
    Ft = ones(T1, 1)
    yt = ones(T1,1)
    statet = zeros(T1,1)

    # initial condition
    i = 100
    Ft[1,] = F[i,]
    yt[1,] = yt[i,]
    statet[1,] = i
    uxt[1,:] = ux[i,:]
    res = zeros(T,2)

    if burn > 0
        draw = [randn(burn); draw]  # Prepend `burn` random values to `draw`
    end
    



end
