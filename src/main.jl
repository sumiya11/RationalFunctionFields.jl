

"""
    global TODOs:
        . ensure no function failes

"""


"""
    Computes several groebner bases of the given G
    while substituting coeffs of the underlying ideal with random points
    Allows to grasp the general structure of the groebner basis
        
    Returns the number of polynomials 
    and the number of field generators (i.e coeficients)
    present in the Groebner basis
"""
function discover_groebner_structure(G::GroebnerEvaluator)
    # generate two substitution points for coeffs of G
    p1, p2 = generate_point(G), generate_point(G)
        
    # compute two groebner bases at that point
    gb1, gb2 = evaluate(G, p1), evaluate(G, p2)
    
    # obtain general info
    structure1 = groebner_structure(gb1)
    structure2 = groebner_structure(gb2)

    if structure1 != structure2
        # :D
        @debug "discovering failed.. recursive call"
        return discover_groebner_structure(G)
    end

    npolys, ncoeffs = structure1
    return npolys, ncoeffs
end


"""
    Returns an array of degrees expected to occure in Groebner basis,
    i-th array entry corresponds to coeff degree with respect to i-th variable

    Gleb: how exactly is array built?
"""
function discover_groebner_degrees(G::GroebnerEvaluator)
    # staring from small amount of interpolation points 
    initial_n = 8
    n = initial_n

    nxs = nvars(G)
    ground = base_ring(G)
    npolys, ncoeffs = discover_groebner_structure(G)
    ncoeffs = sum(ncoeffs)
    
    lightring, = AbstractAlgebra.PolynomialRing(ground, "x")  
    
    # generator for inerpolation sequence
    # sequence_gen^1, sequence_gen^2, sequence_gen^3,...
    sequence_gen = ground(2)
    
    # idx --> max degree
    answer = zeros(Int, nxs)    

    @debug "nvariables " nxs

    for varidx in 1:nxs
        n = initial_n   
        predicted_degrees = [ n for _ in 1:ncoeffs ]
        all_success = false
        
        @debug "handling variable " varidx
        
        while !all_success
            
            # criterion of interpolation success:
            # interpolated degree is far less than the amount of interpolation points
            check_success = (deg) -> 2*(deg + 2) < n

            p = generate_point(G)
            xs = [ deepcopy(p) for _ in 1:n ]
            
            # TODO: sep into a function
            # Gleb: how about broadcast: xs[k] = sequence_gen.^(1:n)  ?
            for k in 1:n
                xs[k][varidx] = rand(ground)
            end
            
            # Gleb: why not @debug ?

            # TODO: do this more effective
            # Gleb: What is this ?!
            ys = map(c -> julia.(c), map(field_generators, map(G, xs)))
            
            all_success = true
            for (j, deg) in enumerate(predicted_degrees)
                xi = [ x[varidx] for x in xs ]
                yi = [ y[j] for y in ys ]
                
                # TODO: fix
                interpolated = lightring(0)

                try
                    interpolated = interpolate_rational_function(lightring, xi, yi)
                catch
                    @warn "rational interpolation failed =("
                    all_success = false
                    break
                end

                predicted_degrees[j] = sum(applytofrac(degree, interpolated))
                
                all_success = all_success && check_success(predicted_degrees[j])
            end

            @debug "for $n points and $varidx variable degrees are " predicted_degrees

            n = n * 2
        end
        
        answer[varidx] = maximum(predicted_degrees)

    end

    answer
end

###############################################################################

"""
    Returns a plain array of coefficients of the given groebner basis
"""
function field_generators(groebner_generators)
    generators = []
    for poly in groebner_generators
        # Gleb: push!(generators, coefficients(poly)...) should work
        for c in coefficients(poly)
            push!(generators, c)
        end
    end
    generators
end

function generators_to_saturated_ideal(genset)
    basepolyring = parent(numerator( first(genset) ))
    nvariables = length(gens(basepolyring))
    ground = base_ring(basepolyring)

    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end

    Fs = map(numerator ∘ (g -> g * Q), genset)

    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = AbstractAlgebra.PolynomialRing(
                                                 basepolyring,
                                                 ystrings,
                                                 ordering=:lex
                                                 )

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    # Gleb: to cancel Fyi and Qy by their gcd
    I = [
        Fyi * Qx - Fxi * Qy
        for (Fyi, Fxi) in zip(Fx, Fy)
    ]

    I, t = saturate(I, Q)

    I, yoverx, basepolyring, nvariables, ground, ystrings, Q, t
end




"""
    Generates the "good" ideal for the given generators set, i.e ideal of form:
    I = <  >
    
    Returns the generated ideal and some auxiliary info
    
"""
# Gleb: As far as I see, this function is always followed by saturation (and this makes perfect sense) 
# I suggest to do the saturation right here
function generators_to_ideal(genset)
    basepolyring = parent(numerator( first(genset) ))
    nvariables = length(gens(basepolyring))
    ground = base_ring(basepolyring)

    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end

    Fs = map(numerator ∘ (g -> g * Q), genset)

    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = AbstractAlgebra.PolynomialRing(
                                                 basepolyring,
                                                 ystrings,
                                                 ordering=:lex
                                                 )

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    # Gleb: to cancel Fyi and Qy by their gcd
    I = [
        Fyi * Qx - Fxi * Qy
        for (Fyi, Fxi) in zip(Fx, Fy)
    ]
    
    # Gleb: how about a named tuple here? Impossible to remember the order
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q
end

"""
    Saturates the given ideal by adding 1 - Q*t to its generators
    where t is the new introduced variable
    
    Returns the saturated ideal and the t variable
"""
function saturate(I, Q)
    R = parent(I[1])
    base = base_ring(R)
    strings = map(String, symbols(R))
    parentring, vs = AbstractAlgebra.PolynomialRing(base, [strings..., "t"], ordering=:lex)

    t = last(vs)
    It = [
        change_parent_ring(f, parentring)[1]
        for f in I
    ]
    sat = 1 - Q*t
    push!(It, sat)
    It, t
end

###############################################################################

"""
    Returns the set of new generators for the field genereted by the given genset
    with no interpolation involved in computing
"""
function naive_new_generating_set(genset)
    
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(genset)

    It, t = saturate(I, Q)
    
    # TODO: ask : why could this naive computation fail for lex ordering 
    #             and work fine for degrevlex
    # Gleb: degrevlex is typically much more efficient, just common wisdom (I do not have a good concise explanation)
    basepolyrings, = Singular.AsEquivalentSingularPolynomialRing(basepolyring)
    yoverxs,  = Singular.PolynomialRing(basepolyrings, [ystrings..., "t"])

    Is = [
          map_coefficients(
              c -> change_base_ring(tosingular(ground), c, parent=basepolyrings),
              f
          )
          for f in It
    ]

    Is = map(F -> change_base_ring(base_ring(yoverxs), F, parent=yoverxs), Is)

    Ideal = Singular.Ideal(yoverxs, Is...)

    gb = collect(gens(Singular.std(Ideal, complete_reduction=true)))

    gb = [
        double_singular2aa(f, base=basepolyring, new_ring=parent(t))
        for f in gb
    ]

    gb_no_t = []
    for f in gb
        if degree(f, t) == 0
            push!(gb_no_t, f)
        end
    end

    gb = [
        change_parent_ring(f, yoverx)
        for f in gb_no_t
    ]
    
    return gb
end

function idealize_and_eval_element(elem, eval_ring, point)
    A, B = numerator(elem), denominator(elem)
    
    A, _ = change_parent_ring(A, eval_ring)
    B, _ = change_parent_ring(B, eval_ring)
    
    println( A )

    Ae = evaluate(A, point)
    Be = evaluate(B, point)

    f = Ae * B - A * Be

    f
end

function idealize_element(elem, basepolyring, basepolyrings, yoverx, yoverxs)
    A, B = numerator(elem), denominator(elem)
    
    ground = base_ring(basepolyring)
    varstrings = string.(symbols(yoverx))

    tmp, = AbstractAlgebra.PolynomialRing(basepolyring, varstrings)

    Ay = change_base_ring(basepolyring, A, parent=tmp)
    By = change_base_ring(basepolyring, B, parent=tmp)
    
    Ay = map_coefficients(c -> c // 1, Ay)
    By = map_coefficients(c -> c // 1, By)
    
    f = A*By - Ay*B
   
    f = map_coefficients(
             c -> change_base_ring(ground, numerator(c), parent=basepolyrings) //
                    change_base_ring(ground, denominator(c), parent=basepolyrings),
             f
    )
    
    f = change_base_ring(
              base_ring(yoverxs),
              f,
              parent=yoverxs)

    f
end

function aa_ideal_to_singular(I)
    basepolyring = parent(first(I))
    nvs = nvars(basepolyring)

    xstrings = ["x$i" for i in 1:nvs]
    ystrings = ["y$i" for i in 1:nvs]
    # Gleb: do you assume that the input variables are always x1, x2, ... ?

    F = base_ring(basepolyring)

    yoverx, yoverxvars = AbstractAlgebra.PolynomialRing(basepolyring, ystrings, ordering=:lex)

    basepolyrings, = Singular.AsEquivalentSingularPolynomialRing(basepolyring)
    fracbases = FractionField(basepolyrings)
    yoverxs, = Singular.PolynomialRing(fracbases,
                             ystrings,
                             ordering=ordering_lp(nvs)*ordering_c())

    Fs = base_ring(basepolyrings)
    
    Is = [
          map_coefficients(c -> change_base_ring(Fs, numerator(c), parent=basepolyrings) // change_base_ring(Fs, denominator(c), parent=basepolyrings), f) for f in I
    ]

    Is = map(c -> change_base_ring(base_ring(yoverxs), c, parent=yoverxs), Is)
    
    Is, yoverx, basepolyring, yoverxs, basepolyrings, fracbases
end




function groebner_structure(gb)
    l = length(gb)
    a = map(f -> length(coefficients(f)), gb)
    l, a
end


function groebner_ideal_to_singular(polys)
    
    yoverx = parent(polys[1])
    varstrings = string.(Singular.symbols(yoverx))
    frac_ring = base_ring(polys[1])
    basepolyring = base_ring(frac_ring)
    ground = base_ring(basepolyring)        

    s_basepolyring, = Singular.AsEquivalentSingularPolynomialRing(basepolyring)
    s_frac = FractionField(s_basepolyring)
    s_yoverx,  = Singular.PolynomialRing(s_frac, varstrings)

    polys = [
        map_coefficients(c -> 
           change_base_ring(ground, numerator(c), parent=s_basepolyring) //
                change_base_ring(ground, denominator(c), parent=s_basepolyring),
           f
        )
        for f in polys
    ]
    
    polys = [
        change_base_ring(base_ring(s_yoverx), f, parent=s_yoverx)
        for f in polys
    ]

    polys, yoverx, basepolyring, s_yoverx, s_basepolyring 
end


"""
    does something

    genset: an array of AbstractAlgebra polynomials over AbstractAlgebra.QQ

    Currently supports 
    AA polynomials + Singular QQ / N_ZpField coefficients
"""
function new_generating_set(genset; modular=true)

    #=
        TODO:
            . ensure evaluated generators are of the same structure
            .
    =#
    
    ### AA polys !!!
    
    modulo = 2^31-1
    FF = Singular.N_ZpField(modulo)
    
    if modular
        genset = modular_reduction(genset, FF)
    end

    # Generating "good" ideal and saturating it
    It, yoverx, basepolyring, nvariables, ground, ystrings, Q, t = generators_to_saturated_ideal(genset)
    
    # Estimating the largest degree of a coeff in the Groebner basis 
    eval_ring, evalvars = Singular.PolynomialRing(
                                        ground,
                                        [ystrings..., "t"],
                          ordering=ordering_lp(nvariables+1)*ordering_c())
    t = evalvars[end]
    
    G = GroebnerEvaluator(It, eval_ring, basepolyring, ground)

    true_structure = discover_groebner_structure(G)
    
    println( "Groebner structure " , true_structure )
    
    exponents = discover_groebner_degrees(G)
    maxexp = maximum(exponents) + 1
    # maxexp = sum(exponents)

    println( "Maximal exponents per variables: $exponents" )
    println( "The largest: $maxexp" )

    # Building a Kronecker substitution

    # why do I even pass npoints to generate_kron...
    npoints = (maxexp + 1)^nvariables + 2
    xs, points = generate_kronecker_points(ground, npoints, maxexp, nvariables)
    println( "The total number of sampling points: $npoints" )
    

    # Evaluating Groebner bases of specializations
    gbs = [
        G(point)
        for point in points
    ]
    
    for gb in gbs
        if groebner_structure(gb) != true_structure
            println(gb)
            @warn "beda"
            @assert false
        end
    end

    gbs_no_t = []
    for gb in gbs
        push!(gbs_no_t, filter(f -> degree(f, t) == 0, gb))
    end
    gbs = gbs_no_t

    
    # Performing univariate interpolation
    uni, = AbstractAlgebra.PolynomialRing(ground, "x")

    ys = []
    for i in 1:npoints
        push!(ys, [Dict() for _ in 1:length(gbs[i])])

        for (j, gg) in enumerate(gbs[i])
            for (ee, c) in zip(exponent_vectors(gg), coefficients(gg))
                ys[i][j][ee] = c
            end
        end

    end
    
    # TODO: reconsider this scary cycle

    answer_e = [Dict() for _ in 1:length(gbs[1])]
    for (j, gg) in enumerate(gbs[1])

        for ev in exponent_vectors(gg)

            yssmall = map(julia, [
                get(y[j], ev, ground(0))
                for y in ys
            ])
            
            f = interpolate_rational_function(
                    uni,
                    xs,
                    yssmall
            )
            
            answer_e[j][ev] = backward_kronecker(f, basepolyring, maxexp)

        end
    end

    # Reconstruct the final Groebner basis and new generators

    gb = []
    R, = AbstractAlgebra.PolynomialRing(
                                FractionField(basepolyring),
                                ystrings,
                                ordering=:lex
                           )

    for etoc in answer_e
        polybuilder = MPolyBuildCtx(R)
        for e in keys(etoc)
            push_term!(polybuilder, etoc[e], e[1:end-1])
        end
        push!(gb, finish(polybuilder))
    end

    if modular
        gb = rational_reconstruction(gb, BigInt(modulo))
    end

    return gb
end
















