

module RationalFunctionFields

    using DataStructures
    using Logging

    using AbstractAlgebra
    import AbstractAlgebra: exponent_vectors, setcoeff!
    import AbstractAlgebra: QQ, divexact, symbols, leading_coefficient, isconstant
    import AbstractAlgebra: coeffs, change_base_ring, gens, base_ring, map_coefficients
    import AbstractAlgebra.Generic
    import AbstractAlgebra.Generic: Frac, MPoly

    import Nemo
    import Nemo: gfp_fmpz_elem, gfp_elem, fmpq_mpoly, gfp_mpoly, fmpz_mpoly
    
    import Singular
    import Singular: std, Ideal
    import Singular: libSingular
    import Singular: spoly, n_Q, ordering_c, ordering_lp
    
    import GroebnerBasis

    import Primes


    include("myeval.jl")
    include("parse.jl")
    include("utils.jl")
    include("interpolation.jl")
    include("modular.jl")
    include("ideals.jl")
    include("groebner.jl")
    # include("inclusion.jl")
    include("algorithms.jl")
    include("RationalFunctionField.jl")

    export RationalFunctionField, simplify!, contains
end




