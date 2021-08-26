
using .RationalFunctionFields: RationalFunctionField, compute_groebner!,
                    contains_using_groebner, contains_randomized

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)


@testset "Basics tests" begin

end

@testset "Contains tests" begin
    FF = Nemo.QQ(2^31 - 1)
    R, (x1, x2) = Nemo.PolynomialRing(FF, ["x1", "x2"])
    # set = [ (a^2 + b^2) // 1, a*b // 1 ]
    # set = [ a*b // (a + b), (a + b) // 1, (a^3 - b) // 1 ]
    set = [
          ( x1^2 + x2^2 ) // R(1),
          ( x1^3 + x2^3 ) // R(1),
          ( x1^4 + x2^4 ) // R(1)
    ]


    FF = RationalFunctionField(set)
    compute_groebner!(FF)  
    
    @test contains_using_groebner( FF, x1*x2 // R(1) ) 
    @test contains_using_groebner( FF, (x1 + x2) // R(1) ) 

    # note that may fail occasionaly
    @test contains_randomized( FF, (x1*x2) // R(1) )
    @test contains_randomized( FF, (x1 + x2 + x1*x2) // R(1) )

end


