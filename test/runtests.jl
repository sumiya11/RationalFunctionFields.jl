using Test
using TestSetExtensions
using Logging

include("../src/RationalFunctionFields.jl")

using .RationalFunctionFields

# Ask me about this
AA = RationalFunctionFields.AbstractAlgebra
Sing = RationalFunctionFields.Singular
Nemo = RationalFunctionFields.Nemo
Primes = RationalFunctionFields.Primes

@info "Testing started"

@testset "All the tests" begin

    @includetests ARGS    

end

@info "All tests OK"
