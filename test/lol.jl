import Singular
import AbstractAlgebra

using .RationalFunctionFields: new_generating_set, generators_to_ideal, GroebnerEvaluator,
                         discover_groebner_structure, discover_groebner_degrees, saturate,
                         naive_new_generating_set, field_generators

logger = Logging.SimpleLogger(stdout, Logging.Debug)
Logging.global_logger(logger)


@testset "lol" begin
    FF = AbstractAlgebra.QQ

    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])

    # true generating coeffs here:
    # -1, x1*x2, -1, -1, x1^2 + x2^2, x1*x2, 1, -x1^2 - x2^2
    set = [
           (x1^2 + x2) // 1,
           (x1*x2) // 1
    ]

    newset = new_generating_set(set)

    println(newset)
    # good test
    @test true
        

end


