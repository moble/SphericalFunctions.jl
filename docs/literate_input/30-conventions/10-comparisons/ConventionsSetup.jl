@testsnippet ConventionsSetup begin

# The generator from which a page draws its samples of angles, created afresh with the same
# seed in every page, so that each page's samples are the same whichever pages have already
# run in the same process
using Random
rng = Random.Xoshiro(1234)

end
