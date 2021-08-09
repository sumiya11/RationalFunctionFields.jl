using Test
using TestSetExtensions


include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields

# Ask me about this
AA = RationalFunctionFields.AbstractAlgebra
Sing = RationalFunctionFields.Singular

@info "Testing started"

@testset "All the tests" begin

    @includetests ARGS    

end

@info "All tests OK"
