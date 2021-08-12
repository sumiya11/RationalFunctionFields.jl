
mutable struct RationalFunctionField
    generating_set       # present always
    groebner_coeffs      # lazy field, not to be addressed directly
    groebner_ideal       # lazy field, not to be addressed directly
end


function RationalFunctionField(genset)
    RationalFunctionField(genset, [], [])
end

function gens(FF::RationalFunctionField)
    return FF.groebner_coeffs
end

function contains_randomized(FF::RationalFunctionField, elem)
    # TODO: merge into a funuction
    I, yoverx, basepolyring, nvariables, ground, ystrings, Q = generators_to_ideal(FF.generating_set)
    It, t = saturate(I, Q)

    @info "" yoverx basepolyring ground

    # Estimating the largest degree of a coeff in the Groebner basis
    eval_ring, evalvars = Singular.PolynomialRing(
                                        ground,
                                        [ystrings..., "t"],
                                        ordering=:lex)
    t = evalvars[end]

    G = GroebnerEvaluator(It, eval_ring)

    p = generate_point(G)
    
    gb = evaluate(G, p)
    x = evaluate(elem, p)
        
    println( gb )
    println( typeof(first(gb)) )
    
    gb = filter(f -> degree(f, t) == 0, gb)

    println( gb )
    println( x, " :: ", typeof(x) )
    
    i = Singular.Ideal(eval_ring, x)
    I = Singular.Ideal(eval_ring, gb)

    return Singular.contains(I, numerator(elem) - x * denominator(elem))
    
end


function contains_using_groebner(FF::RationalFunctionField, elem)
    # elem = A // B
    A, B = numerator(elem), denominator(elem)
    
    Is, yoverx, basepolyring, yoverxs, basepolyrings, fracbases = aa_ideal_to_singular(FF.groebner_ideal)
    
    Ay = change_base_ring(basepolyring, A, parent=yoverx)
    By = change_base_ring(basepolyring, B, parent=yoverx)
    f = change_base_ring(base_ring(yoverxs), A * By - Ay * B, parent=yoverxs)  
    
    fs = Ideal(yoverxs, f)
    Is = Ideal(yoverxs, Is)
    
    println( FF.groebner_coeffs )

    return Singular.contains(Is, fs)
end


function compute_groebner!(
                               FF::RationalFunctionField;
                               backend_algorithm=new_generating_set)
    
    FF.groebner_ideal = backend_algorithm(FF.generating_set)
    FF.groebner_coeffs = field_generators(FF.groebner_ideal)
end


function simple_preprocess(FF::RationalFunctionField)
    generators = FF.groebner_coeffs
    tmp = copy(generators)
    tmp = filter((c) -> !isconstant(c), tmp)
    sort!(tmp, by=(c) -> sum(degrees(c)), rev=true)
    ans = [ first(tmp) ]
    for elem in tmp
        field = RationalFunctionField(ans)
        compute_groebner!(field)
        if !contains(field, elem)
            push!(ans, elem)
        end
    end
    ans
end

function simplify_generators!(
                              FF::RationalFunctionField;
                              backend_algorithm=simple_preprocess)
    FF.groebner_coeffs = backend_algorithm(FF)
end


function contains(FF::RationalFunctionField, elem; proved=true)
    if proved
        contains_using_groebner(FF, elem)
    else
        contains_randomized(FF, elem)
    end
end


###############################################################################

















