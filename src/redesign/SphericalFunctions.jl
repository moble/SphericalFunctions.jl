module Redesign

import Quaternionic
import TestItems: @testitem, @testsnippet
import OffsetArrays


include("WignerDMatrices.jl")
include("WignerMatrix.jl")


function WignerD(R::Quaternionic.Rotor, ℓₘₐₓ::IT; m′ₘₐₓ::IT=ℓₘₐₓ, mₘₐₓ::IT=ℓₘₐₓ) where {IT}
    NT = complex(Quaternionic.basetype(R))
    D = WignerDMatrices(NT, ℓₘₐₓ; m′ₘₐₓ, mₘₐₓ)
    WignerD!(D, R)
end

function WignerD!(D::WignerDMatrices{Complex{FT1}}, R::Quaternionic.Rotor{FT2}) where {FT1, FT2}
    throw(ErrorException("In `WignerD!`, float type of D=$FT1 and of R=$FT2 do not match"))
end

function WignerD!(D::WignerDMatrices{Complex{FT}}, R::Quaternionic.Rotor{FT}) where {FT}
    throw(ErrorException("WignerD! is not yet implemented"))
end

end  # module Redesign
