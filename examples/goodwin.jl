
include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields
using .RationalFunctionFields: load_generators, groebner_ideal,
                            RationalFunctionField, compute_groebner!, gens

Singular = RationalFunctionFields.Singular
Logging = RationalFunctionFields.Logging

logger = Logging.SimpleLogger(stderr, Logging.Debug)
Logging.global_logger(logger)

set = load_generators("../data/Goodwin.txt")

println(set, "\n", length(set))

FF = RationalFunctionField(set)

compute_groebner!(FF)

println( "\nnew gens ", gens(FF) )

