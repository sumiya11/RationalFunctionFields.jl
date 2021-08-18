using .RationalFunctionFields: new_generating_set, generators_to_ideal, GroebnerEvaluator, 
                         discover_groebner_structure, discover_groebner_degrees, saturate,
                         naive_new_generating_set, field_generators, aa_ideal_to_singular,
                         RationalFunctionField, contains_randomized, check_ideal_inclusion

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)


@testset "Preliminary test" begin
    

    R, (x1, x2) = AA.PolynomialRing(Sing.QQ, ["x1", "x2"])

    set = [
           x1 // (x1*x2)
           1 // (x1*x2)
   ]

    # generators_to_ideal test
   I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(set)

    y1, y2 = AA.gens(yoverx)
    x1, x2 = AA.gens(basepolyring)

    @test Q == x1*x2
    # x1*y1*y2 + x1*x2*y1, -y1*y2 + x1*x2
    @test I == [-x1*y1*y2 + x1*x2*y1, -y1*y2 + x1*x2]
    @test ground == Sing.QQ
    @test basepolyring == R
    @test AA.base_ring(yoverx) == basepolyring
    
    # saturate test
    It, t = saturate(I, Q)
    
    f = last(It)
    @test t == last(AA.gens(AA.parent(f)))
    @test f == 1 - x1*x2*t
    
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

    @test Set(coeffs) == Set([x1*x2, -x1*x2, -x1^2 - x2^2, R(1)])

end


@testset "Groebner structure tests" begin

    @info "Groebner structure tests"

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
    eval_ring, = Sing.PolynomialRing(FF, ystrings)
    G = GroebnerEvaluator(I, eval_ring, R, ground)
    npolys, ncoeffs = discover_groebner_structure(G)

    @test npolys == 3 && sum(ncoeffs) == 8

    # over GF
    FF = Sing.N_ZpField(2^31-1)

    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])
    set = [
           (x1^2 + x2^2) // 1,
           (x1*x2) // 1
    ]
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(set)
    eval_ring, = Sing.PolynomialRing(FF, ystrings)
    G = GroebnerEvaluator(I, eval_ring, R, ground)
    
    npolys, ncoeffs = discover_groebner_structure(G)

    @test npolys == 3 && sum(ncoeffs) == 8


    R, (x, y, z) = AA.PolynomialRing(Sing.QQ, ["x", "y", "z"])
    gens = [ (x^2 - z^2) // (x + y), x*y // (x + y), 1//z, z^2 // 1 ] 

    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(gens)
    eval_ring, = Sing.PolynomialRing(Sing.QQ, ystrings)
    G = GroebnerEvaluator(I, eval_ring, R, ground)

    npolys, ncoeffs = discover_groebner_structure(G)
    
    println( npolys, " !!!!!!!!!!!!! ",  ncoeffs )
    
end


@testset "Groebner degree prediction tests" begin
    @info "Groebner degree prediction tests"


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
        eval_ring, = Sing.PolynomialRing(FF, [ystrings..., "t"])
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


function test_groebner_inclusion(set, ismodular)
    newset = new_generating_set(set, modular=ismodular)
    
    @test check_ideal_inclusion(newset, set)

end

function test_field_inclusion(set, ismodular)
    FF = RationalFunctionField(set)    

    for x in set
        @test contains_randomized(FF, x)
    end
end

@testset "Groebner modular and rational are same tests" begin
    FF = Sing.QQ
    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])
    set1 = [
           (x1 + x2) // 1,
           x1*x2 // 1
    ]
    set2 = [
            (x1*x2) // (x1 + 1),
            x2 // 1
    ]

    for set in [set1, set2]
        rational_gb = new_generating_set(set, modular=false)
        finite_gb = new_generating_set(set, modular=true)
        
        @test rational_gb == finite_gb
    end


end

@testset "Groebner main tests" begin
    
    @info "Groebner main tests"

    FF = Sing.QQ
    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])
    set1 = [
           (x1 + x2) // 1,
           x1*x2 // 1
    ]

    set2 = [
            (x1^4 + x2^4) // 1,
            (x1^3 + x2^3) // 1,
            (x1^2 + x2^2) // 1
    ]
    
    for set in [set1, set2]
        for ifmodular in [true, false]
            test_groebner_inclusion(set, ifmodular)
            test_field_inclusion(set, ifmodular)
        end
    end
    
    R, (x1, x2, x3) = AA.PolynomialRing(FF, ["x1", "x2", "x3"])
    # syzygy computation error inside Singular ?
    set3 = [
            (x1 + x2 + x3) // 1,
            (x3 + 2) // (x1 + 1)
    ]
    
    set4 = [
            (x1*x2) // (1),
            1 // (x1 + x2 + x3),
            (2x2 + 3) // 1
    ]
    set5 = [
            (x1 + x2)^2 // 1,
            (x1 + x3)^2 // 1,
            (x2 + x3)^2 // 1
    ]
        
    # we want only modular here since computations over QQ
    # are hard here
    for (i, set) in enumerate([set3, set4, set5])
        for ifmodular in [true]
            test_groebner_inclusion(set, ifmodular)
            test_field_inclusion(set, ifmodular)
        end
    end

end













