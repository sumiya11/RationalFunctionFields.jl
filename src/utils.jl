
"""
    parent_ring_change(poly, new_ring)
Converts a polynomial to a different polynomial ring
Input
  - poly - a polynomial to be converted
  - new_ring - a polynomial ring such that every variable name
      appearing in poly appears among the generators
Output: a polynomial in new_ring "equal" to poly
"""
function parent_ring_change(poly, new_ring)
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any, 1}()
    for u in symbols(old_ring)
        push!(
            var_mapping,
            findfirst(v -> (string(u) == string(v)), symbols(new_ring))
        )
    end
    builder = MPolyBuildCtx(new_ring)
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        new_exp = [0 for _ in gens(new_ring)]
        for i in 1:length(exp)
            if exp[i] != 0
                if var_mapping[i] == nothing
                    throw(Base.ArgumentError("The polynomial contains a variable not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
                end
            end
        end
        push_term!(builder, new_ring.base_ring(coef), new_exp)
    end
    return finish(builder)
end

function parent_ring_change(f::Generic.Frac, new_ring)
    n, d = unpack_fraction(f)
    return parent_ring_change(n, new_ring) // parent_ring_change(d, new_ring)
end


