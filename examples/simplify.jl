
include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields
using .RationalFunctionFields: load_generators, new_generating_set,
             field_generators, compute_groebner!, simplify_generators!,
             RationalFunctionField, field_generators
                                

Singular = RationalFunctionFields.Singular
Logging = RationalFunctionFields.Logging

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)

set = load_generators("../data/LV_gens.txt")

println(set)
FF = RationalFunctionField(set)
compute_groebner!(FF)
simplify_generators!(FF)

println( "\ngenerators ", field_generators(FF) )
