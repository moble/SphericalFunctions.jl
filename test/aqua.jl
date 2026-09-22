@testitem "Aqua quality assurance tests" begin
    using Aqua: Aqua
    # `ambiguities` is on: the package defines methods on `Base.in`, `Base.:*` and friends
    # for its own types, which is exactly where ambiguities creep in.  One did — `in(::Real,
    # ::WignerRange)` against `Base`'s `in(::Integer, ::AbstractUnitRange{<:Integer})` — and
    # went unnoticed while this was `false`.
    Aqua.test_all(SphericalFunctions)
end
