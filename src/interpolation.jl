
using AbstractAlgebra

include("utils.jl")

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


R, (x1, x2) = PolynomialRing(QQ, ["x1", "x2"])
FF = FractionField(R)
RR, (y1, y2) = PolynomialRing(FF, ["y1", "y2"])

# elem
g = ((x1^2 + x2)//(x2+2))*y1 + (x1//x2^2)*y2


# move x2

points = [
    [1, 1],
    [2, 1],
    [3, 1],
]

points = [
    [1, i]
    for i in 1:5
]

# evaluate at different x2
gs = []
for point in points
    gi = map_coefficients(
        t -> evaluate(t, point),
        g)
    push!(gs, gi)
end

UNI, u = PolynomialRing(QQ, "u")

# interpolate
ys = [ coeff(g, 2) for g in gs ]
f = interpolate_rational_function(UNI, [i for i in 1:5], ys)
