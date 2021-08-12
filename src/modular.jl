

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
        return QQ(0, 1)
    end
    if m < 0
        m = -m
    end
    if a < 0
        a = m - a
    end
    if a == 1
        return QQ(1, 1)
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
        return QQ(r, t)
    end

    throw(DomainError(
        :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    ))
end

function rational_reconstruction(f::Nemo.gfp_elem, m)
    rational_reconstruction(BigInt(f.data), m)
end

function rational_reconstruction(f::AbstractAlgebra.GFElem{T}, m) where {T}
    rational_reconstruction(AbstractAlgebra.lift(f), m)
end

function rational_reconstruction(f::Union{MPoly, gfp_mpoly}, m)
    map_coefficients(c -> rational_reconstruction(c, m), f)
end

function rational_reconstruction(f::Frac, m)
    rational_reconstruction(numerator(f), m) //
        rational_reconstruction(denominator(f), m)
end

function rational_reconstruction(A::Array, m)
    map(c -> rational_reconstruction(c, m), A)
end

#------------------------------------------------------------------------------

function modular_reduction(x, field)
    n, d = field(Int(numerator(x))), field(Int(denominator(x)))
    if iszero(d)
        throw(DomainError(
            :($x), "modular reduction of $x (to $field) does not exist"
        ))
    end
    n // d
end

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












