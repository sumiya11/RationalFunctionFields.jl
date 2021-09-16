

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

    # if groebner structures does not match
    if structure1 != structure2
        # :D
        @debug "discovering failed.. recursive call"
        return discover_groebner_structure(G)
    end

    npolys, ncoeffs = structure1
    return npolys, ncoeffs
end


"""
    Computes the Grobner basis of the given ideal at several random points to 
    estimate the variable-wise degrees of the coefficiet polynomials by applying
    univariate interpolation over each variable
    
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
    ncoeffs = sum(ncoeffs)
    
    lightring, = Nemo.PolynomialRing(ground, "x")  
    
    # idx of var --> max degree
    answer = zeros(Int, nxs)    

    for varidx in 1:nxs
        n = initial_n   
        predicted_degrees = [ n for _ in 1:ncoeffs ]
        all_success = false
        
        # @debug "handling variable " varidx
        
        while !all_success
            
            # criterion of interpolation success:
            # interpolated degree is far less
            # than the amount of interpolation points
            check_success = (deg) -> 2*(deg + 2) < n

            p = generate_point(G)
            xs = [ deepcopy(p) for _ in 1:n ]
            
            # TODO: sep into a function
            # Gleb: how about broadcast: xs[k] = sequence_gen.^(1:n)  ?
            # Alex: this won't work..
            for k in 1:n
                xs[k][varidx] = rand(ground)
            end
           
            evaluated = map(G, xs)
            generators = map(obtain_generators, evaluated)
            
            ys = map(c -> ground.(Int.(c)), generators)

            all_success = true
            for (j, deg) in enumerate(predicted_degrees)
                xi = [ x[varidx] for x in xs ]
                yi = [ y[j] for y in ys ]
                
                # TODO: fix
                interpolated = lightring(0)

                try
                    interpolated = interpolate_rational_function(
                                        lightring, xi, yi
                    )
                catch
                    @warn "rational interpolation failed =(\nincrease sample"
                    all_success = false
                    break
                end
                
                num, den = unpackfrac(interpolated)
                predicted_degrees[j] = maximum([degree(num), degree(den)])

                all_success = all_success && check_success(predicted_degrees[j])
            end

            # @debug "for $n points and $varidx variable degrees are " predicted_degrees

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
### AA.coefficients would be better
function obtain_generators(groebner)
    generators = []
    for poly in groebner
        push!(generators, coefficients(poly)...)
    end
    generators
end

"""
    Returns an ideal of polynomials in ys of form
        < Fyi * Qx - Fxi * Qy >
    for each Fi in genset
    Here Q stands for the lcm of the denominators occuring in genset.

    So that generators of the resulting ideal are polynomials,
    but not fractions

    Saturates the ideal with Q by adding 1 - Qt to its generators

"""
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
    yoverx, yoverxvars = Nemo.PolynomialRing(
                                                 basepolyring,
                                                 ystrings,
                                                 ordering=:lex
                                                 )

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    # Gleb: to cancel Fyi and Qy by their gcdi
    # oops
    I = [
        Fyi * Qx - Fxi * Qy
        for (Fyi, Fxi) in zip(Fy, Fx)
    ]
    
    I, t = saturate(I, Q)

    (I=I, yoverx=yoverx, basepolyring=basepolyring, 
     nvariables=nvariables, ground=ground, ystrings=ystrings,
     Q=Q, t=t)
end



# TODO: to delete
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
    yoverx, yoverxvars = Nemo.PolynomialRing(
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
        for (Fyi, Fxi) in zip(Fy, Fx)
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
    parentring, vs = Nemo.PolynomialRing(base, [strings..., "t"], ordering=:lex)

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
function naive_groebner_ideal(genset)
    
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(genset)

    It, t = saturate(I, Q)
    
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
    
    Ae = evaluate(A, point)
    Be = evaluate(B, point)

    f = Ae * B - A * Be

    f
end

function idealize_element(elem, basepolyring, basepolyrings, yoverx, yoverxs)
    A, B = numerator(elem), denominator(elem)
    
    ground = base_ring(basepolyring)
    varstrings = string.(symbols(yoverx))

    tmp, = Nemo.PolynomialRing(basepolyring, varstrings)

    Ay = change_base_ring(basepolyring, A, parent=tmp)
    By = change_base_ring(basepolyring, B, parent=tmp)
    
    Ay = map_coefficients(c -> c // basepolyring(1), Ay)
    By = map_coefficients(c -> c // basepolyring(1), By)
    
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

    yoverx, yoverxvars = Nemo.PolynomialRing(basepolyring, ystrings, ordering=:lex)

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

    s_basepolyring = tosingular(basepolyring)
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
    Computes the Groebner basis of the "good" ideal of the given genset
    interpolating coefficients as rational functions at several numerical points
"""
function groebner_ideal_backend(genset)
    # Generating "good" ideal and saturating it
    
    #=
    ideal_info = generators_to_saturated_ideal(genset)
    It, yoverx, basepolyring, nvariables, ground, ystrings, Q, t = ideal_info.I, ideal_info.yoverx, ideal_info.basepolyring, ideal_info.nvariables, ideal_info.ground, ideal_info.ystrings, ideal_info.Q, ideal_info.t
    =#
    
    I = idealize(genset, saturated=true)
    ground = I.ground    

    singular_ground = tosingular(ground)  
    
    # Estimating the largest degree of a coeff in the Groebner basis
    
    eval_ring, evalvars = Singular.PolynomialRing(
                                        singular_ground,
                                        [I.ystrings..., "t"],
                          ordering=ordering_lp(I.nvariables+1)*ordering_c())
    
    G = GroebnerEvaluator(I.I, eval_ring, I.basepolyring, I.ground)

    true_structure = discover_groebner_structure(G)
    npolys, ncoeffs = true_structure

    @info "Groebner structure: $npolys polys, $(sum(ncoeffs)) coeffs"

    exponents = discover_groebner_degrees(G)
    maxexp = maximum(exponents) + 1

    @info "Maximal exponents per variables: $exponents"
    @info "The largest + 1: $maxexp"
    
    # Building a Kronecker substitution and generating substitution points

    xs, points = generate_good_kronecker_points(ground, exponents, I.nvariables)
    npoints = length(points)

    @info "The total number of sampling points: $npoints"

    @info "Evaluating Groebner bases..."

    # Evaluating Groebner bases of specializations
    gbs = []
    for (i, point) in enumerate(points)
        push!(gbs, evaluate(G, point))
        if i % 1000 == 0
            @info "Computed $i / $npoints Groebner bases.."
            @info last(gbs)
        end
    end
    
    # Ensure bases are of the same shape
    for (gb_i, gb) in enumerate(gbs)
        if groebner_structure(gb) != true_structure
            @warn "Beda: GB of $(groebner_structure(gb)) against $true_structure"
            @warn "Aborting.."
            @assert false
        end
    end

    # Performing univariate interpolation
    @info "Started interpolation.."
    uni, = Nemo.PolynomialRing(ground, "x")

    # Extract coefficients
    ys = [ [ zeros(ground, 0) for _ in 1:ncoeffs[k] ] for k in 1:npolys  ]
    gbcoeffs = [
        [ 
            map(c -> ground(Int(c)), coefficients(gb[i]))
            for i in 1:npolys
        ]
        for gb in gbs
    ]
    for poly_j in 1:npolys
        for coef_i in 1:ncoeffs[poly_j]
            for gb_k in 1:length(gbs)
                push!(
                    ys[poly_j][coef_i],
                    gbcoeffs[gb_k][poly_j][coef_i]
                )
            end
        end
    end

    answer_e = [ Dict() for _ in 1:npolys ]
    for (poly_j, poly) in enumerate(gbs[1])

        @info "Interpolating $poly_j th out of $npolys polynomials.."
        for (c_i, ev) in enumerate(exponent_vectors(poly))

            yscoeff = ys[poly_j][c_i]

            f = interpolate_rational_function(
                    uni,
                    xs,
                    yscoeff
            )
            
            answer_e[poly_j][ev] = backward_good_kronecker(
                                    f, I.basepolyring, exponents
            )
        end
    end

    # Reconstruct the final Groebner basis and new generators

    @info "Reconstructing polynomials.."

    gb = []
    R, = Nemo.PolynomialRing(
                                FractionField(I.basepolyring),
                                I.ystrings,
                                ordering=:lex
                           )

    for etoc in answer_e
        polybuilder = MPolyBuildCtx(R)
        for e in keys(etoc)
            push_term!(polybuilder, etoc[e], e[1:end-1])
        end
        push!(gb, finish(polybuilder))
    end

    gb
end


"""
    Checks if the good ideal composed of the given generators
    contains in the given groebner_basis ideal
"""
function check_ideal_inclusion(groebner_basis, generators)
    Is, yoverx, basepolyring, yoverxs, basepolyrings = groebner_ideal_to_singular(groebner_basis)
    Is = Ideal(yoverxs, Is)

    for f in generators
        f = idealize_element(f, basepolyring, basepolyrings, yoverx, yoverxs)

        fs = Ideal(yoverxs, f)

        if !Singular.contains(Is, fs)
            return false
        end
    end

    return true
end



"""
    Wrapper for the groebner basis computation function

    Performs Modular reduction and Reconstruction, 
    implements CRT to work with big numbers,
    and conducts correctnes checks
    
    Currently supported input
    .    Fractions of Nemo polynomials with Nemo QQ coefficients
"""
function groebner_ideal(
                   initial_genset::Array{T};
                   modular=true,
                   check=false) where {T<:Frac{Nemo.fmpq_mpoly}}

    #=
        TODO:
            . ensure evaluated generators are of the same structure
            .
    =#
    
    # reduction moduli
    moduli = [ 2^30+3 ]

    groebner_bases = []
    
    while true
        
        modulo = last(moduli)
        FF = Nemo.GF(modulo)

        genset = initial_genset
        
        if modular
            @info "Modular reduction.."
            genset = modular_reduction(initial_genset, FF)
        end
        
        gb = groebner_ideal_backend(genset)
        push!(groebner_bases, gb)

        if modular
            @info "CRT reconstruction.."
                
            # let us just hope no coeff vanishes under a small field
            gb = CRT_reconstruction(groebner_bases, BigInt.(moduli))
            
            @info "Rational reconstruction.."
            gb = rational_reconstruction(gb, prod(BigInt.(moduli)))
        
        end
        
        if check
            @info "Running correctness checks.."
            isideal = check_ideal_inclusion(gb, initial_genset)
            @info "Done.\n Is correct ideal: $isideal"

            if isideal ### && isfield ???
                return gb
            end
        
            @warn "Computed *incorrect* \n$gb \nmodulo $modulo"
        
            push!(moduli, Primes.nextprime(BigInt(modulo + 1)))
    
        else
            return gb
        end
    
    end

    @error "how did we end up like this?.."
end


function groebner_ideal(
                   initial_genset::Array{T};
                   modular=true) where {T<:Frac{MPoly{Singular.n_Q}}}
    
    # ?

    strings = string.(gens(parent(numerator(first(initial_genset)))))
    
    subpolyring, = Nemo.PolynomialRing(Nemo.QQ, strings)

    initial_genset = map(
            f -> change_base_ring(Nemo.QQ, numerator(f), parent=subpolyring) //
                change_base_ring(Nemo.QQ, denominator(f), parent=subpolyring),
            initial_genset
    )

    groebner_ideal(initial_genset, modular=modular)
end

































