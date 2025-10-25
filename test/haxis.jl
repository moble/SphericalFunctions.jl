@testitem "HAxis" setup=[EncodeDecode] begin
    using SphericalFunctions.Redesign: HAxis, ℓ, ℓₘᵢₙ, m′ₘᵢₙ, m′ₘₐₓ
    using .EncodeDecode: encode, decode
end
