

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

        if !iszero(Singular.reduce(fs, GroebnerBasis.f4(Is)))
            return false
        end
    end

    return true
end



















