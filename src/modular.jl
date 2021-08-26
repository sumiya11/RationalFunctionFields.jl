

# rational number reconstruction implementation borrowed from CLUE
# and modified a bit to suit the 'Modern Computer Algebra' definitions
# returns a rational r // h of QQ field in a canonical form such that
#   r // h ≡ a (mod m)
#
# let n = max( λ(a), λ(m) ) , where λ(x) is a number of bits for x
# O(n^2)
"""
    rational_reconstruction(a, m)
Rational number reconstruction implementation borrowed from CLUE
and modified a bit to suit the 'Modern Computer Algebra' definitions.
Returns a rational of QQ field in a canonical form that
is congruent a modulo m
a, m are integers
"""
function rational_reconstruction(a::I, m::I) where {I<:Union{Int, BigInt}}
    # restricting to I<:Union{Int, BigInt} for compatibility reasons

    a = mod(a, m)
    if a == 0 || m == 0
        return Nemo.QQ(0, 1)
    end
    if m < 0
        m = -m
    end
    if a < 0
        a = m - a
    end
    if a == 1
        return Nemo.QQ(1, 1)
    end
    bnd = sqrt(float(m) / 2)

    U = I[1, 0, m]
    V = I[0, 1, a]
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end

    t = abs(V[2])
    r = V[3] * sign(V[2])
    # changed from `<= bnd` to `<= m / bnd`
    if t <= m / bnd && gcd(r, t) == 1
        return Nemo.QQ(r, t)
    end

    throw(DomainError(
        :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    ))
end

### Several ad hoc functions for reconstruction

function rational_reconstruction(f::Nemo.gfp_elem, m)
    rational_reconstruction(BigInt(f.data), m)
end

function rational_reconstruction(f::Nemo.fmpz, m)
    rational_reconstruction(BigInt(f), m)
end

function rational_reconstruction(f::AbstractAlgebra.GFElem{T}, m) where {T}
    rational_reconstruction(AbstractAlgebra.lift(f), m)
end

function rational_reconstruction(f::Union{MPoly, gfp_mpoly, fmpz_mpoly}, m)
    map_coefficients(c -> rational_reconstruction(c, m), f)
end

function rational_reconstruction(f::Frac, m)
    rational_reconstruction(numerator(f), m) //
        rational_reconstruction(denominator(f), m)
end

function rational_reconstruction(A::Array, m)
    map(c -> rational_reconstruction(c, m), A)
end

function rational_reconstruction(f::Singular.n_Zp, m)
    rational_reconstruction(BigInt(Int(f)), m)
end

#------------------------------------------------------------------------------

"""
    Coerces the given rational field element into the given finite field
"""
function modular_reduction(x, field)
    n, d = field(Int(numerator(x))), field(Int(denominator(x)))
    if iszero(d)
        throw(DomainError(
            :($x), "modular reduction of $x (to $field) does not exist"
        ))
    end
    n // d
end

### Reduction ad hoc functions

function modular_reduction(f::Union{MPoly, fmpq_mpoly}, field)
    map_coefficients(c -> modular_reduction(c, field), f)
end

function modular_reduction(f::Frac, field)
    modular_reduction(numerator(f), field) //
        modular_reduction(denominator(f), field)
end

function modular_reduction(x::gfp_fmpz_elem, field)
    return field(convert(BigInt, x))
end

function modular_reduction(x::gfp_elem, field)
    return field(x)
end

function modular_reduction(arr::Array, field)
    map(f -> modular_reduction(f, field), arr)
end

#------------------------------------------------------------------------------

"""
    Chinese remainder algorithm
    (Eucleadean one. Ii there a Gauss version of the algorithm btw?..)
"""
function CRT(moduli, remainders)
    """
    The function implements the CR Algorithm mentioned as Algorithm 5.4
    in  'Modern Computer Algebra', 2nd edition,
        Joachim von zur Gathen and Jurgen Gerhard
    """
    
    m = prod(moduli)
    
    ans = BigInt(0)
    for (i, (mi, vi)) in enumerate(zip(moduli, remainders))
        mi_hat = divexact(m, mi)
        _, s, t = Nemo.gcdx(mi_hat, mi)
        ans = BigInt(mod(ans + vi * s * mi_hat, m))
    end
    
    ans
end


"""
    Looks cool
    TODO
"""
function CRT_reconstruction(groebner_bases, moduli)
    onegb = first(groebner_bases)
    onepoly = first(onegb)
    
    basepolyring = base_ring(base_ring(onepoly))
    yoverx = parent(onepoly)
    
    ystring = string.(symbols(yoverx))
    xstring = string.(symbols(basepolyring))
    
    ground = Nemo.ZZ
    inner_ring, = Nemo.PolynomialRing(ground, xstring)
    outer_ring, = Nemo.PolynomialRing(FractionField(inner_ring), ystring)

    npolys = length(onegb)
    ncoeffs = map(length, onegb)    

    new_basis = []
    
    for (polyidx, poly) in zip(1:npolys, onegb)
        outer_f = deepcopy(poly)
        outer_builder = AbstractAlgebra.MPolyBuildCtx(outer_ring)
        
        for (outer_i, outer_e) in zip(1:ncoeffs[polyidx], exponent_vectors(poly))
            inner_frac = coeff(poly, outer_e)
            num_and_denom = [zero(inner_ring), zero(inner_ring)]

            for (jk, locator) in enumerate([numerator, denominator])
                inner_f = locator(inner_frac)
                inner_builder = AbstractAlgebra.MPolyBuildCtx(inner_ring)
            
                for (inner_i, inner_e) in enumerate(exponent_vectors(inner_f))
                    oldcoeffs = [ 
                        coeff(locator(coeff(gb[polyidx], outer_e)), inner_e)
                        for gb in groebner_bases
                    ]
        
                    newcoeff = CRT(moduli, map(x -> BigInt(x.data), oldcoeffs))
                    push_term!(inner_builder, newcoeff, inner_e)
                end

                num_and_denom[jk] = finish(inner_builder)

            end
            
            num, denom = num_and_denom
            push_term!(outer_builder, num // denom, outer_e)
        end    
        
        
        push!(new_basis, finish(outer_builder))
    end
    
    new_basis
end
































