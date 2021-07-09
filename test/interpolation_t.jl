

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









