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



export matchprod
export homeprod
export SurplusVFI



end
