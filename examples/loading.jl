
include("../src/structs.jl")

# load generators set from file 
set = load_generators("../data/verySimple.txt")
set = map(s -> change_base_ring(AbstractAlgebra.QQ, numerator(s)) // change_base_ring(AbstractAlgebra.QQ, denominator(s)), set)

# initialize RationalFunctionField and compute relative groebner basis
FF = RationalFunctionField(set)

compute_groebner!(FF, backend_algorithm=naive_new_generating_set)

# initial generators and a simplified version
println( "generators: ", gens(FF) )

simplify_generators!(FF)

println( "simplified: ", gens(FF) )








