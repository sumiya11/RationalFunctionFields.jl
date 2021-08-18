

"""
    The structure provides convenient computation of 
    a Groebner basis of given ideal at a point.
    
    See evaluate function for more details
"""
mutable struct GroebnerEvaluator
    underlying_ideal
    
    # we want these *not* to be wrapped into Singular
    eval_ring
    coeff_ring
    ground

    saturated
end

"""
    Convenience ctor 1
"""
function GroebnerEvaluator(ideal, eval_ring; saturated=true)
    f = first(ideal)
    GroebnerEvaluator(ideal,
                       eval_ring,
                       base_ring(f),
                       base_ring(base_ring(f)),
                       saturated
    )
end

"""
    Convenience ctor 2
"""
function GroebnerEvaluator(
                           ideal,
                           eval_ring,
                           coeff_ring,
                           ground;
                           saturated=true)

    GroebnerEvaluator(ideal,
                       eval_ring,
                       coeff_ring,
                       ground,
                       saturated
    )
end


"""
    Evaluates all coeffs of the underlying_ideal of G at the given point p,
    then computing the Groebner basis of the resulting ideal
    Returns the computed Groebner basis

    Groebner basis generators are designed to be the elements
    of G.eval_ring over G.ground
    
    Underlying ideal polynomial coefficients live in G.coeff_ring and 
    must agree with substituting the given point p 

    Dependent of the saturating variable "t" polynomials are 
    erased from the resulting Groebner basis
"""
function AbstractAlgebra.evaluate(G::GroebnerEvaluator, p)
    t = last(gens(G.eval_ring))
    
    Is = [
            map_coefficients(c -> evaluate(c, p), f)
            for f in G.underlying_ideal
         ]
    
    Is = [
            change_base_ring(base_ring(G.eval_ring), f, parent=G.eval_ring)
            for f in Is
    ]
    
    ideal = Singular.Ideal(G.eval_ring, Is)
    println(ideal)    
#    gb = GroebnerBasis.f4(ideal, reducegb=1)
    println("why...")
    # this should never happen ideally (but it does!!)
    #if length(gens(gb)) == 1
    if true 
        @warn "F4 failed. Switching to Singular.std"
        gb = Singular.std(ideal, complete_reduction=true)
    end

    gb = collect(gens(gb))
    
    if G.saturated && string(t) == "t"
        gb = filter(f -> degree(f, t) == 0, gb)
    end
    
    #=
    ideal = Singular.Ideal(G.eval_ring, gb)
    gb = collect(gens(Singular.std(ideal, complete_reduction=true)))
    =#

    gb = sort(gb, by=collect ∘ AbstractAlgebra.exponent_vectors)
   
    normalize(f) = !isconstant(f) ? f * inv(leading_coefficient(f)) : f
    gb = map(normalize, gb)    

    gb
end

"""
    Short-cut for evaluate
"""
function (G::GroebnerEvaluator)(p)
    evaluate(G, p)
end


function AbstractAlgebra.nvars(G::GroebnerEvaluator)
    return length(gens(G.coeff_ring))
end

"""
    Returns a random point suitable for evaluation
"""
function generate_point(G::GroebnerEvaluator)
    [ rand(G.ground) for _ in 1:nvars(G) ]
end

function AbstractAlgebra.base_ring(G::GroebnerEvaluator)
    G.ground
end




