
include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields
using .RationalFunctionFields: load_generators, new_generating_set,
                         field_generators

Singular = RationalFunctionFields.Singular
Logging = RationalFunctionFields.Logging

#logger = Logging.SimpleLogger(stderr, Logging.Debug)
#Logging.global_logger(logger)


# load generators set from file 
set = load_generators("../data/Goodwin.txt")

# initialize RationalFunctionField and compute relative groebner basis
println( set )
println( typeof(set[1]) )

gb = new_generating_set(set)
println( "\n", gb )

println( "\n", field_generators(gb) )
