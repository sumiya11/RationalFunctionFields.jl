

mutable struct GroebnerEvaluator
    underlying_ideal
    
    # we want these *not* to be wrapped into Singular
    eval_ring
    coeff_ring
    ground
end

"""
    Convenience ctor
"""
function GroebnerEvaluator(ideal, eval_ring)
    f = first(ideal)
    GroebnerEvaluator(ideal,
                       eval_ring,
                       base_ring(f),
                       base_ring(base_ring(f))
    )
end


"""
    Evaluates all coeffs in the underlying_ideal of G 
    with the given point p, then computing the Groebner basis of the ideal
"""
function AbstractAlgebra.evaluate(G::GroebnerEvaluator, p)
    # TODO: change to GroebnerBasis.jl
    
    Is = [
            map_coefficients(c -> evaluate(c, p), f)
            for f in G.underlying_ideal
         ]
    
    # @info "" first(Is) parent(first(Is)) base_ring(first(Is))
    # @info "" G.ground G.eval_ring G.coeff_ring
    
    Is = [
            change_base_ring(base_ring(G.eval_ring), f, parent=G.eval_ring)
            for f in Is
    ]

    ideal = Singular.Ideal(G.eval_ring, Is)
    gb = collect(gens(Singular.std(ideal, complete_reduction=true)))
    
    gb
end

# shortcut for evaluate 
function (G::GroebnerEvaluator)(p)
    evaluate(G, p)
end


function AbstractAlgebra.nvars(G::GroebnerEvaluator)
    return length(gens(G.coeff_ring))
end

function generate_point(G::GroebnerEvaluator)
    [ rand(G.ground) for _ in 1:nvars(G) ]
end

function AbstractAlgebra.base_ring(G::GroebnerEvaluator)
    return G.ground
end










