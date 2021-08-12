using .RationalFunctionFields: new_generating_set, generators_to_ideal, GroebnerEvaluator, 
                         discover_groebner_structure, discover_groebner_degrees, saturate,
                         naive_new_generating_set, field_generators, aa_ideal_to_singular 

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)


@testset "Preliminary test" begin
    

    R, (x1, x2) = AA.PolynomialRing(Sing.QQ, ["x1", "x2"])

    set = [
           x1*(x1*x2) // (x1*x2)^2
           1 // (x1*x2)
   ]

    # generators_to_ideal test
   I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(set)

    y1, y2 = AA.gens(yoverx)
    x1, x2 = AA.gens(basepolyring)

    @test Q == x1*x2
    @test I == [-y1^2*y2 + x1^2*x2, -y1*y2 + x1*x2]
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

    # over GF
    FF = Sing.GF(2^31-1)

    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])
    set = [
           (x1^2 + x2^2) // 1,
           (x1*x2) // 1
    ]
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(set)
    eval_ring, = Singular.PolynomialRing(FF, ystrings)
    G = GroebnerEvaluator(I, eval_ring, R, ground)
    
    npolys, ncoeffs = discover_groebner_structure(G)

    @test npolys == 3 && ncoeffs == 8


    ##### hmmm
    # Gleb: ??

    R, (x, y, z) = AA.PolynomialRing(Sing.QQ, ["x", "y", "z"])
    gens = [ (x^2 - z^2) // (x + y), x*y // (x + y), 1//z, z^2 // 1 ] 

    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(gens)
    eval_ring, = Singular.PolynomialRing(Sing.QQ, ystrings)
    G = GroebnerEvaluator(I, eval_ring, R, ground)

    npolys, ncoeffs = discover_groebner_structure(G)
    
    println( npolys, " !!!!!!!!!!!!! ",  ncoeffs )
    
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

#=
@testset "Groebner small tests" begin
    
    # AA polys, Singular.QQ
    
    FF = Singular.QQ
    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])

    set = [
           (x1 + x2) // 1,
           (x1*x2) // 1
    ]
    newset = new_generating_set(set, modular=false)

    println(newset)


    # AA polys, Singular.QQ + modular
    
    newset = new_generating_set(set)

    println(newset) 

end
=#

@testset "Groebner main tests" begin
    
    # AA polys, Singular.QQ

    FF = Singular.QQ
    R, (x1, x2) = AA.PolynomialRing(FF, ["x1", "x2"])
    set = [
           (x1 + x2) // 1,
           x1*x2 // 1
    ]

    newset = new_generating_set(set)
    
    println( AA.parent(newset[1]) )
    println( AA.base_ring(newset[1]) )
    println( AA.base_ring(AA.base_ring(newset[1])) )
    println( typeof(AA.base_ring(AA.base_ring(AA.base_ring(newset[1]))) ))


    Is, _, _, yoverxs, basepolyrings, _ =  aa_ideal_to_singular(newset)
    
    Rs, (xs1, xs2) = Sing.PolynomialRing(FF, ["x1", "x2"])
    # Gleb: why is this again here?
    sets = [
        (xs1 + xs2) // 1,
        xs1*xs2 // 1
    ]
    i, yoverx2, basepolyring2, nvariables2, ground2, ystrings2, Q2 = generators_to_ideal(sets)
    
    Fs = AA.base_ring(basepolyrings)
    is = [
          AA.map_coefficients(c -> AA.change_base_ring(Fs, c, parent=basepolyrings), f) for f in i
    ]
    is = map(c -> AA.change_base_ring(AA.base_ring(yoverxs), c, parent=yoverxs), is)


    I = Sing.Ideal(yoverxs, Is)
    i = Sing.Ideal(yoverxs, is)
    
    
    println(I)
    println(i)

    println(Sing.contains(I, i))
end











