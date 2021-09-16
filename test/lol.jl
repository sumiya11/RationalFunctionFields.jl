
using .RationalFunctionFields

logger = Logging.SimpleLogger(stdout, Logging.Debug)
Logging.global_logger(logger)


@testset "lol" begin
    FF = Nemo.QQ

    R, (x1, x2) = Nemo.PolynomialRing(FF, ["x1", "x2"])

    set = [
          ( x1^2 + x2^2 ) // R(1),
          ( x1 ) // R(1),
    ]

    FF = RationalFunctionField(set)
    
    println( RationalFunctionFields.contains(FF, x1 // R(1), p=1) )
    
    println( RationalFunctionFields.contains(FF, x2^2 // R(1), p=0.99) )


    println( RationalFunctionFields.groebner_ideal(set, check=true) )
end







