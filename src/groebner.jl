


mutable struct groebner_evaluator
    underlying_ideal
    eval_ring
    ground
end


"""
    Evaluates all coeffs in the underlying_ideal of G 
    with the given point p, then computing the Groebner basis of the ideal
"""
function AbstractAlgebra.evaluate(G::groebner_evaluator, p)
    # TODO: change to GroebnerBasis.jl
    
    Is = [
            map_coefficients(c -> evaluate(c, p), f)
            for f in G.underlying_ideal
         ]

    Is = [
            change_base_ring(tosingular(G.ground), f, parent=G.eval_ring)
            for f in Is
    ]

    ideal = Singular.Ideal(G.eval_ring, Is)
    gb = Singular.std(ideal, complete_reduction=true)

    gb
end













