module Robin2011_RepPackage

using Distributions
using Copulas
using LinearAlgebra
using Distributions

"""
matchprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686)

Takes the worker and aggregate state grids, as well as the parameters B and C, and returns the matrix of match productivities.

Default parameter values are those used in Robin (2011).
"""

function matchprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686)
    return y * (B * x .+ C)'
end

"""
homeprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686, α::Number = 0.64, z0::Number = 0.766752794650811)

Takes the worker and aggregate state grids, as well as the parameters B, C, α, and returns the matrix of home production.

Default parameter values are those used in Robin (2011).
"""

function homeprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686, α::Number = 0.64, z0::Number = 0.766752794650811)
    return z0 .+ α * (matchprod(x, y; B, C) .- z0)
end

"""
SurplusVFI(p::Matrix, z::Matrix, Π::Matrix; β::Number = 0.9466)

Performs Value Function Iteration on the Surplus Function, returns the resulting Surplus function matrix.

Takes as inputs the matrix of match productivities p, the matrix of home production z, the markov transition matrix Π, and the parameter β as a keyword argument.

β should be set to equal (1-δ)/(1+r).
"""
function SurplusVFI(p::Matrix, z::Matrix, Π::Matrix; β::Number = 0.946603693905558)
    S = (I(length(Π[1,:])) - β .* Π )\(p - z)
    e = norm(S - max.(S, 0), 2)

    while e > 0.0001
        S1 = S
        S = p - z + β * Π * max.(S, 0)
        e = norm(S - S1, 2)
    end

    return S
end

include("wages.jl")

include("dynamics.jl")

include("grids.jl")

include("estCrit.jl")

include("params_in.jl")

export params_estimated
export estCrit
export wage_dens_path
export SurplusVFI
export WageVFI
export matchprod
export homeprod
export unemp_path

end
