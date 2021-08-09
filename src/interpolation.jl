
function lagrange_interpolant(R, xs, i)
    n = length(xs)
    x = AbstractAlgebra.gen(R)

    f = prod(
        j == i ? 1 : (x - xs[j])
        for j in 1:n
    )
    c = prod(
        j == i ? 1 : (xs[i] - xs[j])
        for j in 1:n
    )
    
    AbstractAlgebra.map_coefficients(t -> t // c, f)
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

    n = AbstractAlgebra.degree(m)
    k = div(n, 2)
    FF = AbstractAlgebra.parent(g)

    U = (FF(1), FF(0), m)
    V = (FF(0), FF(1), g)
    while AbstractAlgebra.degree(V[3]) >= k
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end

    r, t = V[3], V[2]


    if gcd(r, t) == FF(1) && AbstractAlgebra.degree(t) <= n - k
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

    x = AbstractAlgebra.gen(R)

    g = interpolate_polynomial(R, xs, ys)
    m = prod(x - xs[i] for i in 1:length(xs))
    polynomial_reconstruction(g, m)
end


function kronecker_substitution(f)
    r, t = numerator(f), denominator(f)
   

end


function interpolate_multivariate_rational_function(MR, R, xs, ys)
   
    f = interpolate_rational_function(R, xs, ys)
    
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

function random_linear_shift(ground_ring, n)
    functors = []
    inverses = []
    
    for i in 1:n
        a, b = ground_ring(rand(1:10000)), ground_ring(rand(1:10000))
        f  = ( x -> x + b ) 
        fi = ( x -> x - b )
        
        push!(functors, f)
        push!(inverses, fi)
    end
    
    return functors, inverses
end


"""
    Given a callable functor f tries to find degrees of univariate polynomials
    A, B approximating f as a rational function f = A // B
    
    Numeric computations are performed in the given ground field
"""
function predict_degrees(f, ground)
    # staring from small amount of points..
    n = 8

    # = deg(A), deg(B) initially
    predicted_degrees = (n, n)

    polyring, = AbstractAlgebra.PolynomialRing(ground, "u")    
    
    # "safety condition" to assume f is interpolated correctly indeed
    # failes on first iteration
    while 4*(sum(predicted_degrees) + 2) > n  
        xs = [ rand(ground) for _ in 1:n  ]
        ys = f.(xs)
            
        interpolated = interpolate_rational_function(polyring, xs, ys)    
        
        predicted_degrees = map(degree, (numerator(interpolated), denominator(interpolated)))   
        
        @debug "for $n points degrees are " predicted_degrees

        n = n * 2
    end
    
    predicted_degrees
end





















