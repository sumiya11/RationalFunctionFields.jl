using Test
using TestSetExtensions


include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields



@info "Testing started"

@testset "All the tests" begin

    @includetests [:utils_tests, :interpolation_tests, :groebner_tests]

end

@info "All tests OK"
