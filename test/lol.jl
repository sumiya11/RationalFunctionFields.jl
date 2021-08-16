
using .RationalFunctionFields: new_generating_set, generators_to_ideal, GroebnerEvaluator,
                         discover_groebner_structure, discover_groebner_degrees, saturate,
                         naive_new_generating_set, field_generators

logger = Logging.SimpleLogger(stdout, Logging.Debug)
Logging.global_logger(logger)


@testset "lol" begin
    FF = Sing.QQ

    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])

    set = [
          ( x1^2 + x2^2 ) // 1,
          ( x1^3 + x2^3 ) // 1,
          ( x1^4 + x2^4 ) // 1
    ]

    newset = new_generating_set(set)

    println(newset)
    println(typeof(newset[1]))
    println("##############")

    newset = new_generating_set(set, modular=false)

    println(newset)
    println(typeof(newset[1]))
    println(field_generators(newset))

end







