

import AbstractAlgebra
import Singular

import Singular.libSingular
import Singular: spoly, n_Q
import AbstractAlgebra: coeffs, change_base_ring, gens, base_ring, map_coefficients
import AbstractAlgebra.Generic
import AbstractAlgebra.Generic: Frac, MPoly
import Nemo


function unknown2known(u)
    libSingular.julia(libSingular.cast_number_to_void(u.ptr))    
end

function singular2aa(poly::Singular.spoly{T}; base=false, new_ring=false) where {T}
    nvariables = length(gens(parent(poly)))
    xstrings = ["x$i" for i in 1:nvariables]
    if base == false
        base = base_ring(poly)
    end
    if new_ring == false
        new_ring, = AbstractAlgebra.PolynomialRing(base, xstrings)
    end
    change_base_ring(base, poly, parent=new_ring)
end

function double_singular2aa(poly::spoly{Singular.n_unknown{spoly{T}}}; base=false, new_ring=false) where {T}
    outer_change = singular2aa(poly)

    basebase = base_ring(parent(unknown2known(collect(coeffs(poly))[1])))

    nvariables = length(gens(parent(outer_change)))
    ystrings = ["y$i" for i in 1:nvariables]
    new_ring, = AbstractAlgebra.PolynomialRing(basebase, ystrings)

    inner_change = map_coefficients(
                      c -> singular2aa(unknown2known(c), base=basebase, new_ring=base),
                      outer_change)
    inner_change
end


function aa2singular(poly::MPoly{T}; base=false, new_ring=false) where {T}
    nvariables = length(gens(parent(poly)))
    xstrings = ["x$i" for i in 1:nvariables]
    if base == false
        base = base_ring(poly)
    end
    if new_ring == false
        new_ring, = Singular.PolynomialRing(base, xstrings)
    end
    change_base_ring(base, poly, parent=new_ring)    
end

function double_aa2singular(poly::MPoly{T}; base=false, new_ring=false) where {T}
    nvariables = length(gens(parent(poly)))
    ystrings = ["y$i" for i in 1:nvariables]
    xstrings = ["x$i" for i in 1:nvariables]
    basebase = base_ring(base_ring(parent(poly)))
    if base == false
        base, = Singular.PolynomialRing(basebase, xstrings)
    end
    if new_ring == false
        new_ring, = Singular.PolynomialRing(base, ystrings)
    end
    change_base_ring(base_ring(new_ring), poly, parent=new_ring)
end

function Nemo.degree(f::AbstractAlgebra.Generic.Frac{T}) where {T}
    return max(degree(denominator(f)), degree(numerator(f)))
end

function Nemo.isconstant(f::AbstractAlgebra.Generic.Frac{T}) where {T}
    return isconstant(denominator(f)) && isconstant(numerator(f))
end

function Nemo.degrees(f::AbstractAlgebra.Generic.Frac{T}) where {T}
    return degrees(denominator(f)) + degrees(numerator(f))
end

function Base.length(x::Frac{T}) where {T}
    return max(length(numerator(x)), length(denominator(x)))
end

function tosingular(F::AbstractAlgebra.Generic.Rationals)
    Singular.QQ
end

function tosingular(F)
    F
end

function toaa(::Singular.Rationals)
    AbstractAlgebra.Generic.QQ
end

function iota(n)
    return [i for i in 1:n]
end

###############################################################################




