

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

    eval_ring, evalvars = Singular.PolynomialRing(
                              ground,
                              [ystrings..., "t"],
                              # Gleb: why and what is this? To discuss
                              ordering=ordering_lp(nvariables)*ordering_c())

    G = GroebnerEvaluator(It, eval_ring)

    p = generate_point(G)

    gb = evaluate(G, p)
    elem = idealize_and_eval_element(elem, eval_ring, p)
    
    @debug "Evaluated eLement is $elem"
    @debug "Evaluated groebner basis is $gb"

<<<<<<< HEAD
    I = Singular.Ideal(eval_ring, gb)

=======
    @info "Groebner basis is $gb"
    I = Singular.Ideal(eval_ring, gb)
>>>>>>> 6fd74be1ac5fce6df9d128e2b7417cb0ac3ed90a
    return iszero(Singular.reduce(elem, GroebnerBasis.f4(I)))
end


function contains_using_groebner(FF::RationalFunctionField, elem)

    Is, yoverx, basepolyring, yoverxs, basepolyrings = groebner_ideal_to_singular(FF.groebner_ideal)

    f = idealize_element(elem, basepolyring, basepolyrings, yoverx, yoverxs)

    fs = Ideal(yoverxs, f)
    Is = Ideal(yoverxs, Is)
    
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
























