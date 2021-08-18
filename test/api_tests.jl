
using .RationalFunctionFields: RationalFunctionField, compute_groebner!,
                    contains_using_groebner, contains_randomized

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)


@testset "Basics tests" begin

end

@testset "Contains tests" begin
    FF = Sing.QQ
    R, (a, b) = AA.PolynomialRing(FF, ["a", "b"])
    set = [ (a^2 + b^2) // 1, a*b // 1 ]

    FF = RationalFunctionField(set)
    compute_groebner!(FF)  
    
    @test contains_using_groebner( FF, a*b // 1 ) 
    @test contains_using_groebner( FF, (a^2 + a*b + b^2) // 1 ) 
        
    # note that may fail occasionaly
    @test contains_randomized( FF, a*b // 1 )
    @test contains_randomized( FF, (a^2 + b^2)^3 // 1 )

end


