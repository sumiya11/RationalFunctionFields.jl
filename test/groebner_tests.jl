
import Singular
import AbstractAlgebra

using .RationalFunctionFields: new_generating_set, generators_to_ideal, GroebnerEvaluator, 
                         discover_groebner_structure, discover_groebner_degrees, saturate,
                         naive_new_generating_set, field_generators 

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)


@testset "Preliminaries test" begin

end

# essentially a sanity check
@testset "Naive generation tests" begin

    R, (x1, x2) = AA.PolynomialRing(Sing.QQ, ["x1", "x2"])

    set = [
           (x1^2 + x2^2) // 1,
           (x1*x2) // 1
   ]

    gb = naive_new_generating_set(set)
    coeffs = field_generators(gb)

    @test Set(coeffs) == Set([R(-1), x1^2 + x2^2, x1*x2, R(1), -x1^2 - x2^2])

end


@testset "Groebner structure tests" begin

    # ogo, chto-to novoe
    FF = Sing.QQ

    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])

    # true generating coeffs here:
    # -1, x1*x2, -1, -1, x1^2 + x2^2, x1*x2, 1, -x1^2 - x2^2
    set = [
           (x1^2 + x2^2) // 1,
           (x1*x2) // 1
    ]

    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(set)

    eval_ring, = Singular.PolynomialRing(FF, ystrings)

    G = GroebnerEvaluator(I, eval_ring, R, ground)

    npolys, ncoeffs = discover_groebner_structure(G)

    @test npolys == 3 && ncoeffs == 8

end


@testset "Groebner degree prediction tests" begin
    # now to the degree prediction tests

    FF = Sing.QQ

    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])

    set1 = [
        (x1^2 + x2^2) // 1,
        (x1*x2) // 1
    ]

    for (i, set) in enumerate([ set1 ])
        I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(set)
        It, t = saturate(I, Q)
        eval_ring, = Singular.PolynomialRing(FF, [ystrings..., "t"])
        coeff_ring = basepolyring

        G = GroebnerEvaluator(It, eval_ring, coeff_ring, ground)

        predicted_degrees = discover_groebner_degrees(G)

        gb = naive_new_generating_set(set)
        coeffs = field_generators(gb)
        true_degrees = [
            maximum([ AA.degree(f, var) for f in coeffs ])
            for var in AA.gens(coeff_ring)
        ]

        @test predicted_degrees == true_degrees

    end

end


@testset "Groebner main algorithm tests" begin
    
    # AA polys, Singular.QQ
    
    FF = Singular.QQ
    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])

    set = [
           (x1^2 + x2^2) // 1,
           (x1*x2) // 1
    ]
    newset = new_generating_set(set, modular=false)

    println(newset)


    # AA polys, Singular.QQ + modular
    
    newset = new_generating_set(set)

    println(newset) 

end













