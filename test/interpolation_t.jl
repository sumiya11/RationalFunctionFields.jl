

include("../src/interpolation.jl")


function test_rational_function_interpolation(f)
    
    R = parent(denominator(f))
    
    # hope we have multivariate ring here
    vs = collect(gens(R))
    nvs = length(vs)

    npoints = 100

    varidx = 1

    points = [ [ one(QQ) for _ in 1:nvs ] for _ in 1:npoints ]
    for i in 1:npoints
        points[i][varidx] = QQ(i)
    end

    xs = []
    ys = []
    for i in 1:npoints
        push!(xs, points[i][varidx])
        push!(ys, evaluate(numerator(f), points[i]) // evaluate(denominator(f), points[i]))
    end

    uniring, x = PolynomialRing(QQ, "x")

    inter = interpolate_rational_function(uniring, xs, ys)
        
    vs = [x1, x2]

    degf = degree(numerator(f), vs[varidx]), degree(denominator(f), vs[varidx])
    degi = degree(numerator(inter)), degree(denominator(inter))

    println(degi, " : ", degf)

    @assert degi == degf
end


R, (x1, x2, x3) = PolynomialRing(QQ, ["x1", "x2", "x3"])

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


###############################################################################


R, (x1, x2, x3) = PolynomialRing(QQ, ["x1", "x2", "x3"])

# f = (x1^2*x2^2 + x1*x2 + x1) // (x1 + 1234*x2 + 1)
f = (x1^2*x2 + x3) // (x2^2 + x3^2)

npoints = 150
ndeg = 3
nvs = length(gens(R))

points = [ [ BigInt(j)^i for i in [ (ndeg+1)^k for k in 1:nvs ] ] for j in 1:npoints ]

xs = []
ys = []

shifts, inverses = random_linear_shift(QQ, length(points[1]))
fp = [ fj(pj) for (fj, pj) in zip(shifts, gens(R)) ]
g = evaluate(numerator(f), fp) // evaluate(denominator(f), fp)

# println(f, " --> ")
# println(g)

for i in 1:npoints
    p = points[i]

    push!(xs, p[1])
    push!(ys,
       evaluate(numerator(g), p) // evaluate(denominator(g), p)
    )
end

uniring, x = PolynomialRing(QQ, "x")

# println(xs, "\n", ys, "\n-----")

inter = interpolate_multivariate_rational_function(R, uniring, xs, ys, 3)

# println("interpolated = ", inter)

lol = [R(0), R(0)]
for (i, D) in [ [1, numerator(inter)], [2, denominator(inter)] ]
    polybuilder = MPolyBuildCtx(R)
    for (e, c) in zip(0:degree(D), coefficients(D))
        if c == 0
            continue
        end
        push_term!(polybuilder, c, decompose_by_degrees(e, ndeg+1, nvs))
    end
    lol[i] = finish(polybuilder)
end

F = lol[1] // lol[2]

# println("back kronecker F = ", F)

shiftedvariables = [ fj(pj) for (fj, pj) in zip(inverses, gens(R)) ]

F = evaluate(numerator(F), shiftedvariables) // evaluate(denominator(F), shiftedvariables)

F = leading_coefficient(numerator(F)) *inv(leading_coefficient(denominator(F))) * (numerator(F) * inv(leading_coefficient(numerator(F)))) // (denominator(F) * inv(leading_coefficient(denominator(F))))

# println("back linear F = ", F)

@assert abs(F // f) == 1

























