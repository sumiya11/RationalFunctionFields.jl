

function lagrange_interpolant(R, xs, i)
    n = length(xs)
    x = AbstractAlgebra.gen(R)

    f = R(prod(
        j == i ? 1 : (x - xs[j])
        for j in 1:n
    ))
    
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


"""
    For given polynomial g and m computes the rational fraction r//t,
    s.t. r//t approximates g modulo m 
    
        r == t*g (mod m)
    
    and
        deg(r) + deg(t) <= deg(m)

    Throws if ??? TODO
"""
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


"""
    Interpolates rational function in the given ring R with
    the given interpolation points xs and ys by
    interpolating as a polynomial first
    and reconstructing into a polynomial fraction after
"""
function interpolate_rational_function(R, xs, ys)
    # The function implements algorithm derived from Corollary 5.20 of
    #   "Modern Computer Algebra", second edition,
    #   Joachim von zur Gathen, Jürgen Gerhard

    x = AbstractAlgebra.gen(R)
    g = Nemo.interpolate(R, xs, ys)

    # TODO: this computation can be done faster
    m = prod(x - xs[i] for i in 1:length(xs))
    
    polynomial_reconstruction(g, m)
end


# TODO: why is it still here?
function interpolate_multivariate_rational_function(MR, R, xs, ys)
   
    f = interpolate_rational_function(R, xs, ys)
    
    f   
end

"""
    Returns reversed vector of base-adic representation of n
"""
function decompose_by_degrees(n, base, len)
    powers = [base^(i-1) for i in 1:len]
    exps = zeros(Int, len)
    for i in 1:len
        exps[len - i + 1] = div(n, powers[len - i + 1])
        n = n - exps[len - i + 1] * powers[len - i + 1]
    end
    exps
end


# TODO : not used
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
"""
# TODO: not used
function predict_degrees(f)
    # staring from small amount of points..
    n = 8

    # = deg(A), deg(B) initially
    predicted_degrees = (n, n)

    polyring, = AbstractAlgebra.PolynomialRing(ground, "u")    
    
    # "safety condition" to assume f is interpolated correctly indeed
    # failes on first iteration
    while 4*(sum(predicted_degrees) + 2) > n  
        xs = [ rand(f) for _ in 1:n  ]
        ys = f.(xs)
            
        interpolated = interpolate_rational_function(polyring, xs, ys)    
        
        predicted_degrees = map(degree, (numerator(interpolated), denominator(interpolated)))   
        
        @debug "for $n points degrees are " predicted_degrees

        n = n * 2
    end
    
    predicted_degrees
end



"""
    Same as the function above, key difference is that here
    f returns an array of results on each evaluation
"""
# TODO: not used
function simultaneous_predict_degrees(f)
    # staring from small amount of points..
    n = 8
    
    polyring = base_ring(f)
    nfrom, nto = arity(f)

    predicted_degrees = [ (n, n) for _ in 1:nto ]
    already_interpolated = 0    

    while already_interpolated < nto 
        
        4*(sum(values(predicted_degrees)) + 2) > n
        xs = [ rand(f) for _ in 1:n  ]
        ys = f.(xs)

        for (j, deg) in enumerate(predicted_degrees)
            if 4*(sum(deg) + 2) < n
                continue
            end
            yparticular = [y[j] for y in ys]
            interpolated = interpolate_rational_function(polyring, xs, yparticular)
            
            predicted_degrees[j] = map(degree, (numerator(interpolated), denominator(interpolated)))

        end

        @debug "for $n points degrees are " predicted_degrees

        n = n * 2
    end

    predicted_degrees
end


"""
    
"""
# TODO: not used
function backward_kronecker(
                     f,
                     target_ring,
                     maxexp) where {T}

    polybuilder = MPolyBuildCtx(target_ring)
    nvariables = length(gens(target_ring)) 
    

    # TODO: fix 0:degree(f)
    for (e, c) in zip(0:degree(f), coefficients(f))
        if iszero(c)
            continue
        end
        push_term!(
            polybuilder,
            c,
            decompose_by_degrees(e, maxexp, nvariables)
        )
    end

    finish(polybuilder)
end


function decompose_by_good_degrees(e, prods, n)
    ans = zeros(Int, n)

    for i in 1:n
        ans[n - i + 1] = div(e, prods[n - i + 1])
        e = mod(e, prods[n - i + 1])
    end
    
    ans
end

function backward_good_kronecker(
                     f,
                     target_ring,
                     exps) where {T}
    
    prods = deepcopy(exps)
    prods[1] = 1
    for i in 2:length(prods)
        prods[i] = prods[i - 1] * (exps[i - 1] + 1)
    end

    polybuilder = MPolyBuildCtx(target_ring)
    nvariables = length(gens(target_ring))

    for (e, c) in zip(0:degree(f), coefficients(f))
        if iszero(c)
            continue
        end
        push_term!(
            polybuilder,
            c,
            decompose_by_good_degrees(e, prods, nvariables)
        )
    end

    ans = finish(polybuilder)
    ans
end


function backward_kronecker(
                  f::AbstractAlgebra.Generic.Frac{T},
                  target_ring,
                  maxexp) where {T}

    backward_kronecker(numerator(f), target_ring, maxexp) //
        backward_kronecker(denominator(f), target_ring, maxexp)
end

function backward_good_kronecker(
                  f::AbstractAlgebra.Generic.Frac{T},
                  target_ring,
                  maxexp) where {T}

    backward_good_kronecker(numerator(f), target_ring, maxexp) //
        backward_good_kronecker(denominator(f), target_ring, maxexp)
end

function generate_kronecker_points(ground, maxexp, nvariables)
    npoints = (maxexp^nvariables + 2) * 2
    xs = [rand(ground) for _ in 1:npoints]
    points = [
        [ ground(j)^i for i in [ (maxexp)^k for k in 0:(nvariables - 1) ] ]
          for j in xs
    ]
    xs, points
end

function generate_good_kronecker_points(ground, exps, nvariables)
    npoints = prod( map(x -> x + 1, exps) ) * 2
   
    prods = deepcopy(exps)
    # univariate case? =(
    prods[1] = 1
    for i in 2:length(prods)
        prods[i] = prods[i - 1] * (exps[i - 1] + 1)
    end
    
    xs = Set([rand(ground) for _ in 1:npoints])
    while length(xs) < npoints
        xs = union!(xs, Set([rand(ground) for _ in 1:(npoints-length(xs))]))
    end
    xs = collect(xs)

    points = [
        [ p^prods[k] for k in 1:nvariables ]
        for p in xs
    ]
    
    xs, points
end



















