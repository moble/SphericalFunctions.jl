@testitem "Aqua quality assurance tests" begin
    using Aqua: Aqua
    Aqua.test_all(
        SphericalFunctions;
        ambiguities=false,
        #stale_deps=(;ignore=[:Requires])  # Need Requires on Julia <1.9; not loaded on ≥1.9
    )
end
