

"""
    The structure provides convenient computation of 
    a Groebner basis of given ideal at a point.
    
    See evaluate function for more details
"""

### Probably we want to dispatch on coefficients type
mutable struct GroebnerEvaluator
    ideal::IdealContext
end

"""
    Convenience ctor 1
    
    ?? convenience
"""
function GroebnerEvaluator(ideal)
    GroebnerEvaluator(ideal)
end

#=
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
=#


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
    context = G.ideal
    
    I = context.I
    ground = context.ground
    evalring = context.evalring    

    singular_ground = base_ring(evalring)

    # TODO: change
    lift = typeof(p[1]) <: Nemo.gfp_elem ? x -> x.data : x -> x
        
    Is = [
            map_coefficients(c -> singular_ground(lift(evaluate(c, p))), f)
            for f in I
         ]

    Is = [
            change_base_ring(base_ring(evalring), f, parent=evalring)
            for f in Is
    ]
    
    ideal = Singular.Ideal(evalring, Is)
    # @info "I am gpoing to compute GB of $ideal"
    gb = GroebnerBasis.f4(ideal, reducegb=0, monorder=:lex)
    # this should never happen ideally (but it does!!)
    # @info gb
    if length(gens(gb)) == 1
        @info gb
        @warn "F4 failed. Switching to Singular.std"
        gb = Singular.std(ideal, complete_reduction=true)
    end

    gb = collect(gens(gb))
    
    # ideal = Singular.Ideal(eval_ring, gb)
    # gb = collect(gens(Singular.std(ideal, complete_reduction=true)))

    # @info "After reduction $gb"

    if context.saturated
        t = last(gens(evalring))
        gb = filter(f -> degree(f, t) == 0, gb)
    end
    
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
    return length(gens(G.ideal.basepolyring))
end

"""
    Returns a random point suitable for evaluation
"""
function generate_point(G::GroebnerEvaluator; M=Inf)
    if M == Inf
        [ rand(G.ideal.ground) for _ in 1:nvars(G) ]
    else
        [ G.ideal.ground(rand(0:M)) for _ in 1:nvars(G) ]
    end
end

function AbstractAlgebra.base_ring(G::GroebnerEvaluator)
    G.ideal.ground
end




