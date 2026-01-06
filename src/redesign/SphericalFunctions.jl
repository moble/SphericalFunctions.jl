module Redesign

import Quaternionic
import TestItems: @testitem, @testsnippet
import FixedSizeArrays: FixedSizeArrayDefault, FixedSizeVectorDefault,  FixedSizeArray,
    FixedSizeVector

# TODO: remove once this moves to SphericalFunctions proper
import SphericalFunctions: ComplexPowers


include("Wigner/wigner.jl")


# function WignerD(R::Quaternionic.Rotor, ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ) where {IT}
#     NT = complex(Quaternionic.basetype(R))
#     D = WignerDMatrices(NT, ℓₘₐₓ, m′ₘₐₓ)
#     WignerD!(D, R)
# end

# function WignerD!(D::WignerDMatrices{Complex{FT1}}, R::Quaternionic.Rotor{FT2}) where {FT1, FT2}
#     R1 = Quaternionic.Rotor{FT1}(R)
#     WignerD!(D, R1)
# end

# function WignerD!(D::WignerDMatrices{Complex{FT}}, R::Quaternionic.Rotor{FT}) where {FT}
#     error("WignerD! is not yet implemented")
# end

end  # module Redesign
