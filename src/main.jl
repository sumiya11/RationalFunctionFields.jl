

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
    npolys1, npolys2 = length(gb1), length(gb2)
    ncoeffs1, ncoeffs2 = map(length ∘ field_generators, [gb1, gb2])
    
    if npolys1 != npolys2 || ncoeffs1 != ncoeffs2
        # :D
        @debug "discovering failed.. recursive call"
        return discover_groebner_structure(G)
    end

    npolys1, ncoeffs1
end


"""
    Returns an array of degrees expected to occure in Groebner basis,
    i-th array entry corresponds to coeff degree with respect to i-th variable
"""
function discover_groebner_degrees(G::GroebnerEvaluator)
    # staring from small amount of interpolation points 
    initial_n = 8
    n = initial_n

    nxs = nvars(G)
    ground = base_ring(G)
    npolys, ncoeffs = discover_groebner_structure(G)
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
            for k in 1:n
                xs[k][varidx] = sequence_gen^k
            end

            # TODO: do this more effective
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
        for c in coefficients(poly)
            push!(generators, c)
        end
    end
    generators
end

"""
    Generates the "good" ideal for the given generators set, i.e ideal of form:
    I = <  >
    
    Returns the generated ideal and some auxiliary info
    
"""
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

    I = [
        Fyi*Qx - Fxi*Qy
        for (Fyi, Fxi) in zip(Fx, Fy)
    ]
    
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

# TODO: to be deleted
function exponents_new_generating_set(genset)

    I, basepolyring, nvariables, ground, ystrings = generators_to_ideal(genset)

    # substitute + groebner n times

    ybasering,  = AbstractAlgebra.PolynomialRing(ground, ystrings)

    # TODO: restructure so that interpolation is attampted after each GB computation
    n = 30

    ans = Dict()

    for (varidx, var) in enumerate(gens(basepolyring))

        gbs = []

        # generating interpolation points while fixing all variables but var

        # TODO: take other points random
        # TODO ?: custom random for each field 
        rand_sample = [ ground( rand(-1000:1000) ) for i in 1:nvariables ] 
        points = [
            deepcopy(rand_sample)
            for _ in 1:n
        ]
        for i in 1:length(points)
            points[i][varidx] = ground(i)
        end
        
        lightring, = Singular.PolynomialRing(tosingular(ground), ystrings)

        for point in points
            gb = evaluate_gb_at_point(I, point, lightring, ground)
            push!(gbs, gb)
        end

        gbs = [
            [
                change_base_ring(ground, f, parent=basepolyring)
                for f in gens(gb)
            ]
            for gb in gbs
        ]

        ints = []

        for (j, g) in enumerate(gbs[1])
            push!(ints, Dict())

            lcoeff = leading_coefficient(g)
            g = map_coefficients(c -> c // lcoeff, g)

            exps = collect(exponent_vectors(g))
            for (k, c) in enumerate(monomials(g))

                xs = [ point[varidx] for point in points ]
                ys = [ coeff(gb[j], k) for gb in gbs ]

                UNI, = AbstractAlgebra.PolynomialRing(ground, "x")
                f = interpolate_rational_function(UNI, xs, ys)

                ints[j][exps[k]] = f
            end
        end

        ans[var] = [ ints,  gbs[1] ]

    end

    return ans
end

"""
    Returns the set of new generators for the field genereted by the given genset
    with no interpolation involved in computing
"""
function naive_new_generating_set(genset)
    
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(genset)

    It, t = saturate(I, Q)
    
    # TODO: ask : why could this naive computation fail for lex ordering 
    #             and work fine for degrevlex
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


# TODO: to be deleted
function compare_degrees(true_basis, interpolated_basis)

    total = 0
    hit = 0
    
    
    for (varidx, var) in enumerate(keys(interpolated_basis))

        cffs, basis = interpolated_basis[var]

        true_coeffs = []


        for f in gens(true_basis)
            push!(true_coeffs, Dict())


            f = double_singular2aa(f)
            lcoeff = leading_coefficient(f)
            f = map_coefficients(c -> c // lcoeff, f)

            true_exps = [collect(exponent_vectors(m)) for m in monomials(f)]

            for (i, c) in enumerate(coeffs(f))
                true_coeffs[end][true_exps[i]...] = c
            end
        end

        Rxs = base_ring(true_coeffs[1][[0, 0]])
        xsvars = collect(gens(Rxs))
        
        ground = base_ring(Rxs)

        xstrings = ["x$i" for i in 1:length(xsvars)]
        Rxsaa, myxs = AbstractAlgebra.PolynomialRing(ground, xstrings, ordering=:lex)

        x = xsvars[varidx]
        myx = myxs[varidx]


        myeverything = Dict{Any, Set}()
        trueverything = Dict{Any, Set}()

        for (f, g) in zip(cffs, true_coeffs)
            for e in keys(f)
                if !haskey(myeverything, e)
                    myeverything[e] = Set()
                end
                if !haskey(trueverything, e)
                    trueverything[e] = Set()
                end

                myterm = f[e]

                trueterm = 1
                if !haskey(g, e)
                    for gg in true_coeffs
                        if haskey(gg, e)
                            trueterm = gg[e]
                            break
                        end
                    end
                else
                    trueterm = g[e]
                end

                mydeg = max(degree(numerator(myterm)),
                           degree(denominator(myterm)))


                xs = collect(gens(parent(numerator(trueterm))))

                truedeg = max(degree(numerator(trueterm), xs[varidx]),
                           degree(denominator(trueterm), xs[varidx]))

                push!(myeverything[e], mydeg)
                push!(trueverything[e], truedeg)
            end
        end


        for (key, val) in myeverything
            for dg in val
                if haskey(trueverything, key)
                    if dg in trueverything[key]
                        hit += 1
                    else
                        @warn "bedabedaogorchenie"
                        # return hit, total
                    end
                else
                    @warn "beda"
                    return hit, total
                end
                total += 1
            end
        end
    end


    return hit, total
end

# TODO: to be deleted
function display_generating_exponents(genset)
    inter = exponents_new_generating_set(genset)

    answer = [ Dict() for _ in 1:length(first(values(inter))[1]) ]

    for v in keys(inter)
        exps, basis = inter[v]

        for (idx, g) in enumerate(exps)
            for ee in keys(g)
                if !haskey(answer[idx], ee)
                    answer[idx][ee] = 0
                end
                dg = degree(g[ee])
                answer[idx][ee] += dg
            end
        end
    end

    return answer
end


function aa_ideal_to_singular(I)
    basepolyring = parent(first(I))
    nvs = nvars(basepolyring)

    xstrings = ["x$i" for i in 1:nvs]
    ystrings = ["y$i" for i in 1:nvs]

    F = base_ring(basepolyring)

    yoverx, yoverxvars = AbstractAlgebra.PolynomialRing(basepolyring, ystrings)

    basepolyrings, = Singular.AsEquivalentSingularPolynomialRing(basepolyring)
    fracbases = FractionField(basepolyrings)
    yoverxs, = Singular.PolynomialRing(fracbases, ystrings)

    Fs = base_ring(basepolyrings)
    
    Is = [
          map_coefficients(c -> change_base_ring(Fs, numerator(c), parent=basepolyrings) // change_base_ring(Fs, denominator(c), parent=basepolyrings), f) for f in I
    ]

    Is = map(c -> change_base_ring(base_ring(yoverxs), c, parent=yoverxs), Is)
    
    Is, yoverx, basepolyring, yoverxs, basepolyrings, fracbases
end


"""
    does something

    genset: an array of AbstractAlgebra polynomials over AbstractAlgebra.QQ
"""
function new_generating_set(genset; modular=true)
    
    #=
        TODO:
            . ensure evaluated generators are of the same structure
            .
    =#
    
    ### AA polys !!!
    
    modulo = 2^31 - 1
    FF = Singular.GF(modulo)
    
    if modular
        genset = modular_reduction(genset, FF)
    end

    # Generating "good" ideal and saturating it
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(genset)
    It, t = saturate(I, Q)
    
    @info "" yoverx basepolyring ground

    # Estimating the largest degree of a coeff in the Groebner basis 
    eval_ring, evalvars = Singular.PolynomialRing(
                                        ground,
                                        [ystrings..., "t"],
                                        ordering=:lex)
    t = evalvars[end]
    
    G = GroebnerEvaluator(It, eval_ring)

    exponents = discover_groebner_degrees(G)
    maxexp = sum(exponents)

    @debug "Maximal exponents per variables: $exponents"
    @debug "The largest: $maxexp"

    # Building a Kronecker substitution

    npoints = (maxexp + 1)^nvariables + 2
    xs, points = generate_kronecker_points(ground, npoints, maxexp, nvariables)
    println( points )
    @debug "The total number of sampling points: $npoints"
    

    # Evaluating Groebner bases of specializations
    gbs = [
        G(point)
        for point in points
    ]
    
    extended_ring, vvs = AbstractAlgebra.PolynomialRing(
                                                 ground,
                                                 [ystrings..., "t"],
                                                 ordering=:lex
                                         )

    gbs_no_t = []
    for gb in gbs
        push!(gbs_no_t, filter(f -> degree(f, t) == 0, gb))
    end
    gbs = gbs_no_t


    # Performing univariate interpolation

    # 
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

    answer_e = [Dict() for _ in 1:length(gbs[1])]
    for (j, gg) in enumerate(gbs[1])

        for ev in exponent_vectors(gg)

            yssmall = map(julia, [
                haskey(y[j], ev) ? y[j][ev] : ground(0)
                for y in ys
            ])
            

            f = interpolate_multivariate_rational_function(
                    basepolyring,
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
















