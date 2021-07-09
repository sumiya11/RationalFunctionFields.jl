
using AbstractAlgebra

include("utils.jl")

GROUND = GF(2^31 - 1)


function lagrange_interpolant(R, xs, i)
    n = length(xs)
    x = gen(R)

    f = prod(
        j == i ? 1 : (x - xs[j])
        for j in 1:n
    )
    c = prod(
        j == i ? 1 : (xs[i] - xs[j])
        for j in 1:n
    )
    
    map_coefficients(t -> t // c, f)
end

# univariate
# returns :: AbstractAlgebra.Generic.Poly
function interpolate_polynomial(R, xs, ys)
    f = sum(
        ys[i] * lagrange_interpolant(R, xs, i)
        for i in 1:length(xs)
    )
end


# Returns the minimal polynomial of the sequence `h` modulo `m`,
# the sequence `h` is formed as the given polynomial coefficients
# If n is max(deg(h), deg(m))
# O(n^2)
function polynomial_reconstruction(g, m)
    # The function implements Algorithm 12.9 from
    #   "Modern Computer Algebra", second edition,
    #   Joachim von zur Gathen, Jürgen Gerhard

    n = degree(m)
    k = div(n, 2)
    FF = parent(g)

    U = (FF(1), FF(0), m)
    V = (FF(0), FF(1), g)
    while degree(V[3]) >= k
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end

    r, t = V[3], V[2]


    if gcd(r, t) == FF(1) && degree(t) <= n - k
        return r // t
    end

    throw(DomainError(
        g, "reconstruction of $g (mod $m) does not exist"
    ))
end

# univariate
# returns :: AbstractAlgebra.Generic.MPoly
function interpolate_rational_function(R, xs, ys)
    # The function implements algorithm derived from Corollary 5.20 of
    #   "Modern Computer Algebra", second edition,
    #   Joachim von zur Gathen, Jürgen Gerhard

    x = gen(R)

    g = interpolate_polynomial(R, xs, ys)
    m = prod(x - xs[i] for i in 1:length(xs))
    polynomial_reconstruction(g, m)
end


function kronecker_substitution(f)
    r, t = numerator(f), denominator(f)
   

end


function interpolate_multivariate_rational_function(MR, R, xs, ys, termdeg)
   
    f = interpolate_rational_function(R, xs, ys)
    
    vars = collect(gens(MR))

#    for t in terms(f)
        
 #   end

    println(f)
    
    f   
end

function decompose_by_degrees(n, base, len)
    powers = [base^(i-1) for i in 1:len]
    exps = zeros(Int, len)
    for i in 1:len
        exps[len - i + 1] = div(n, powers[len - i + 1])
        n = n - exps[len - i + 1] * powers[len - i + 1]
    end
    exps
end

function random_linear_shift(n)
    functors = []
    inverses = []
    
    for i in 1:n
        a, b = GROUND(rand(1:100)), GROUND(rand(1:100))
        f  = ( x -> x + b ) 
        fi = ( x -> x - b )
        
        push!(functors, f)
        push!(inverses, fi)
    end
    
    return functors, inverses
end

#=
R, (x1, x2, x3) = PolynomialRing(GROUND, ["x1", "x2", "x3"])

# f = (x1^2*x2^2 + x1*x2 + x1) // (x1 + 1234*x2 + 1)
f = (x1^2*x2 + x3) // (x2^2 + x3^2)

npoints = 150
ndeg = 3
nvs = length(gens(R))

points = [ [ BigInt(j)^i for i in [ (ndeg+1)^k for k in 1:nvs ] ] for j in 1:npoints ]

xs = []
ys = []

shifts, inverses = random_linear_shift(length(points[1]))
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

uniring, x = PolynomialRing(GROUND, "x")

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




=#














