

import AbstractAlgebra

include("interpolation.jl")

# just for convenience
AA = AbstractAlgebra


mutable struct polynomial_interpolation
    groebner
    polyring
    xs
    ys
end

function eat_point(P::polynomial_interpolation, p)
    
    x, y = p
    
    push!(P.xs, x)
    push!(P.ys, y)
        
    if length(P.xs) == 40
        f = interpolate_polynomial(P.polyring, P.xs, P.ys)
        return f
    end
    
    return P.polyring(0)
end

function update(P::polynomial_interpolation, f)
    P.xs = []
    P.ys = []
end


mutable struct rational_interpolation
    polynomial
    predicted_degrees
    init
    power
    current
end

function generate_point(R::rational_interpolation)
    
    base = R.polynomial.groebner.base
    R.current = R.current * R.init
    R.power += 1
    
    dA, dB = R.predicted_degrees
    
    return R.current .* [base(2), base(3)] .* (base(5))^(dA - dB)
end

function interpolate_coeff(R::rational_interpolation)
    
    # power -> interpolated coeff
    polyring = R.polynomial.polyring
    
    A, B = polyring(0), polyring(0)
    
    while iszero(A) && iszero(B)
        nevals = sum(R.predicted_degrees) + 2
        xs = [ generate_point(R) for i in 1:nevals ]
        ys = [ evaluate(R.polynomial.groebner, xs[i]) for i in 1:nevals ]
        xs = [ x[1] for x in xs ]
        Q = interpolate_rational_function(R.polynomial.polyring, xs, ys)
        
        Ai, Bi = map(leading_coefficient, [numerator(Q), denominator(Q)])
        println(Q)
        
        A = eat_point(R.polynomial, Ai)
        B = eat_point(R.polynomial, Bi)
    end
    
    return A // B
end

function update(R::rational_interpolation, f)
    R.power = 0
    R.current = R.init^R.power
end

mutable struct groebner_evaluation
    f
    index
    base
end

function update(G::groebner_evaluation, f)
    G.index += 1
end

function evaluate(G::groebner_evaluation, p)
    return G.f(p...)[G.index]
end


FF = AA.GF(2^31 - 1)

# emulates groebner basis coeffs computation
f = (x, y) -> FF.([ (x^2 + y) // (x + 11), (x*y + 1) // (8x^3 + 13) ])


function predict_degrees(f)
    #return [(2, 1), (2, 3)]
    return (2, 1)
end

function interpolate(f)
    # predict shape (i.e positions of coeffs) of f
    predicted_degrees = predict_degrees(f)
    ncoeffs = 2  
    
    base = FF

    interpolation_ring, = AA.PolynomialRing(base, "x")
    
    groebner = groebner_evaluation(f, 1, base)
    
    polynomial = polynomial_interpolation(groebner, interpolation_ring, [], [])
    
    rational = rational_interpolation(polynomial, predicted_degrees, base(2), 0, base(2)^0)
    
    i = 0
    while i < ncoeffs
        c = interpolate_coeff(rational)
        @info "interpolated " c
        i += 1
    end
    
end


interpolate(f)










