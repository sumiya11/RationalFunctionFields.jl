
using .RationalFunctionFields: load_generators 

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)


@testset "Loading from file tests" begin

    set = load_generators("../data/Simple.txt")
    
    R = parent(first(set))
    x, y, z = AA.gens(AA.base_ring(R))

    @test set[1] == (x^2 - z^2) // (x + y)
    @test set[2] == (x*y) // (x + y)
    @test set[3] == (1) // z
    @test set[4] == z^3 // z


end

