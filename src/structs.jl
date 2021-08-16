

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
    # TODO: use Theorem 3.3 from the draft

    # TODO: merge into a funuction
    It, yoverx, basepolyring, nvariables, ground, ystrings, Q, t = generators_to_saturated_ideal(FF.generating_set)

    # Estimating the largest degree of a coeff in the Groebner basis
    eval_ring, evalvars = Singular.PolynomialRing(
                              ground,
                              [ystrings..., "t"],
                              ordering=ordering_lp(nvariables)*ordering_c())
    t = evalvars[end]

    G = GroebnerEvaluator(It, eval_ring)

    p = generate_point(G)

    gb = evaluate(G, p)
    elem = idealize_and_eval_element(elem, eval_ring, p)

    println( gb )
    println( typeof(first(gb)) )
    println( elem , typeof(elem) )
    println(parent(elem) == parent(gb[1]))
            
    i = Singular.Ideal(eval_ring, elem)
    I = Singular.Ideal(eval_ring, gb)

    println(Singular.std(I))

    return Singular.contains(I, i)
end


function contains_using_groebner(FF::RationalFunctionField, elem)

    Is, yoverx, basepolyring, yoverxs, basepolyrings = groebner_ideal_to_singular(FF.groebner_ideal)

    f = idealize_element(elem, basepolyring, basepolyrings, yoverx, yoverxs)

    fs = Ideal(yoverxs, f)
    Is = Ideal(yoverxs, Is)
    
    println(Is, " | ", fs)
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
























