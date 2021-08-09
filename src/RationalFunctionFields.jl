

module RationalFunctionFields

    using DataStructures
    using Logging

    using AbstractAlgebra
    import AbstractAlgebra: QQ, divexact, symbols
    import AbstractAlgebra: coeffs, change_base_ring, gens, base_ring, map_coefficients
    import AbstractAlgebra.Generic
    import AbstractAlgebra.Generic: Frac, MPoly

    import Nemo
   
    import Singular
    import Singular: std, Ideal
    import Singular: libSingular
    import Singular: spoly, n_Q

    include("myeval.jl")
    include("parse.jl")
    # include("modular.jl")
    include("utils.jl")
    include("interpolation.jl")
    include("groebner.jl")
    include("main.jl")
    include("structs.jl")

end
