
include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields
using .RationalFunctionFields: load_generators, new_generating_set,
                         field_generators

Singular = RationalFunctionFields.Singular
Logging = RationalFunctionFields.Logging

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)

set = load_generators("../data/Sontag_gen.txt")

gb = new_generating_set(set)

println( "\ngb ", gb )
println( "\ngenerators ", field_generators(gb) )
