
using .RationalFunctionFields: new_generating_set, generators_to_ideal, GroebnerEvaluator,
                         discover_groebner_structure, discover_groebner_degrees, saturate,
                         naive_new_generating_set, field_generators

logger = Logging.SimpleLogger(stdout, Logging.Debug)
Logging.global_logger(logger)


@testset "lol" begin
    FF = Nemo.QQ

    R, (x1, x2) = Nemo.PolynomialRing(FF, ["x1", "x2"])

    set = [
          ( x1^2 + x2^2 ) // R(1),
          ( x1^3 + x2^3 ) // R(1),
          ( x1^4 + x2^4 ) // R(1)
    ]

    newset1 = new_generating_set(set)
    # newset2 = new_generating_set(set, modular=false)
    
    println(newset1)
    # println(newset2)
    println("#######")
    println(field_generators(newset1))
    # println(field_generators(newset2))

    
    
end







