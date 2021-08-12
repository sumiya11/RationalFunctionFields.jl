
include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields
using .RationalFunctionFields: load_generators, new_generating_set,
                         field_generators

Singular = RationalFunctionFields.Singular
Logging = RationalFunctionFields.Logging

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)


# load generators set from file 
set = load_generators("../data/Goodwin.txt")

set = map(s -> Singular.change_base_ring(Singular.QQ, numerator(s)) // Singular.change_base_ring(Singular.QQ, denominator(s)), set)

# initialize RationalFunctionField and compute relative groebner basis
println( set )

gb = new_generating_set(set)
println( "\n", gb )

println( "\n", field_generators(gb) )






