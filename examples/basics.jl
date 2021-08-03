
include("../src/structs.jl")


# create polyring and a generators set
R, (a, b) = AbstractAlgebra.PolynomialRing(Singular.QQ, ["a", "b"])
set = [ (a + b) // 1, b // 1 ]

# initialize RationalFunctionField and compute relative groebner basis
FF = RationalFunctionField(set)
compute_groebner!(FF)  # default method is used

# initial generators and a simplified version
println( "generators: ", gens(FF) )

simplify_generators!(FF)

println( "simplified: ", gens(FF) )

# compute_groebner!(FF, backend_algorithm=naive_new_generating_set)






