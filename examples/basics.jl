include("../src/RationalFunctionFields.jl")

using AbstractAlgebra
using Logging
using Singular

using .RationalFunctionFields: simplify_generators!, compute_groebner!, RationalFunctionField

logger = SimpleLogger(stdout, Logging.Debug)

# create polyring and a generators set
R, (a, b) = AbstractAlgebra.PolynomialRing(Singular.QQ, ["a", "b"])
set = [ (a^2 + b^2) // 1, (a^3 + b^3) // 1, (a^4 + b^4) // 1 ]

# initialize RationalFunctionField and compute relative groebner basis
FF = RationalFunctionField(set)
compute_groebner!(FF)  # default method is used

# initial generators and a simplified version
println( "generators: ", gens(FF) )

simplify_generators!(FF)

println( "simplified: ", gens(FF) )

# compute_groebner!(FF, backend_algorithm=naive_new_generating_set)






