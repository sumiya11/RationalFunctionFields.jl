
if !isdefined(Main, :testset)
    using Test
    using TestSetExtensions
end

include("../src/RationalFunctionFields.jl")
using .RationalFunctionFields: interpolate_rational_function, interpolate_multivariate_rational_function, random_linear_shift, decompose_by_degrees, predict_degrees


import Singular
import AbstractAlgebra


AA = AbstractAlgebra

function test_rational_function_interpolation(f)
    
    R = parent(denominator(f))
    
    # hope we have multivariate ring here
    vs = collect(AA.gens(R))
    nvs = length(vs)

    npoints = 100

    varidx = 1

    points = [ [ one(AA.QQ) for _ in 1:nvs ] for _ in 1:npoints ]
    for i in 1:npoints
        points[i][varidx] = AA.QQ(i)
    end

    xs = []
    ys = []
    for i in 1:npoints
        push!(xs, points[i][varidx])
        push!(ys, AA.evaluate(numerator(f), points[i]) // AA.evaluate(denominator(f), points[i]))
    end

    uniring, x = AA.PolynomialRing(AA.QQ, "x")

    inter = interpolate_rational_function(uniring, xs, ys)
        
    degf = AA.degree(numerator(f), vs[varidx]), AA.degree(denominator(f), vs[varidx])
    degi = AA.degree(numerator(inter)), AA.degree(denominator(inter))

    @test degi == degf
end

@testset "Rational Interpolation tests" begin

    R, (x1, x2, x3) = AA.PolynomialRing(AA.QQ, ["x1", "x2", "x3"])

    fs = [
      (x1 + x2) // x2,
      x1^5 // (x2 + x1 + x3),
      (x2 + x3) // (x1^3 + x1*x2 + x1*x3),
      (x1^4*x2^3 + x1^4*x3 + x2) // x1,
      (x1^6 + x2) // (x1^9 + x3),
      (x1^36 + x3^12 + x1^23 - x2^2 - x2*x3^2 + x1^2*x2^3+(2*x3^4)) // (x1^2 + x2 + x1*x2 + x3)
    ]

    for f in fs
        test_rational_function_interpolation(f)
    end

end

###############################################################################


@testset "Interpolation + shift tests" begin

    R, (x1, x2, x3) = AA.PolynomialRing(AA.QQ, ["x1", "x2", "x3"])

    f = (x1^2*x2 + x3) // (x2 + 3)

    npoints = 100
    ndeg = 3
    nvs = length(AA.gens(R))

    points = [ [ BigInt(j)^i for i in [ (ndeg+1)^k for k in 1:nvs ] ] for j in 1:npoints ]

    xs = []
    ys = []

    shifts, inverses = random_linear_shift(AA.QQ, length(points[1]))
    fp = [ fj(pj) for (fj, pj) in zip(shifts, AA.gens(R)) ]
    g = AA.evaluate(numerator(f), fp) // AA.evaluate(denominator(f), fp)


    for i in 1:npoints
        p = points[i]

        push!(xs, p[1])
        push!(ys,
            AA.evaluate(numerator(g), p) // AA.evaluate(denominator(g), p)
        )
    end

    uniring, x = AA.PolynomialRing(AA.QQ, "x")

    inter = interpolate_multivariate_rational_function(R, uniring, xs, ys)

    lol = [R(0), R(0)]
    for (i, D) in [ [1, numerator(inter)], [2, denominator(inter)] ]
        polybuilder = AA.MPolyBuildCtx(R)
        for (e, c) in zip(0:AA.degree(D), AA.coefficients(D))
            if c == 0
                continue
            end
            AA.push_term!(polybuilder, c, decompose_by_degrees(e, ndeg+1, nvs))
        end
        lol[i] = AA.finish(polybuilder)
    end

    F = lol[1] // lol[2]

    shiftedvariables = [ fj(pj) for (fj, pj) in zip(inverses, AA.gens(R)) ]

    F = AA.evaluate(numerator(F), shiftedvariables) // AA.evaluate(denominator(F), shiftedvariables)

    F = AA.leading_coefficient(numerator(F)) *inv(AA.leading_coefficient(denominator(F))) * (numerator(F) * inv(AA.leading_coefficient(numerator(F)))) // (denominator(F) * inv(AA.leading_coefficient(denominator(F))))


    @test isone(denominator(F // f)) && isone(numerator(F // f))

end


@testset "Degree prediction tests" begin
    
    
    fs0 = [
        (x) -> 2,
        (x) -> x,
        (x) -> 1 // x,
        (x) -> (x^2 + 11) // (x + 111)
    ]
    answers0 = [
        (0, 0),
        (1, 0),
        (0, 1),
        (2, 1)
    ]
        
    fs = [
        (x) -> (x^3 + 1),
        (x) -> (x^16 - 23x^11 + 123),
        (x) -> (x^2 // (x + 1)^4),
        (x) -> ((x - 1) // (x^18 + x^17 - 123)),
        (x) -> ((5x^113 + 2x^55 + 3x^33 + 1) // (x^11 + x + 22))
    ]


    answers = [ 
        (3, 0),
        (16, 0),
        (2, 4),
        (1, 18),
        (113, 11)
    ]


    FF = AA.GF(2^31-1)

    for (i, f) in enumerate(fs0)
        d = predict_degrees(f, AA.QQ)
        @test answers0[i] == d
    end
    
    for (i, f) in enumerate(fs)
        d = predict_degrees(f, FF)
        @test answers[i] == d
    end


end


#=

# for future 

using Logging

logger = Logging.SimpleLogger(stderr, Logging.Debug)

Logging.global_logger(logger)
@debug "lol"


=#


















