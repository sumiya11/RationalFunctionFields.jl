

function unpackfrac(q)
    numerator(q), denominator(q)
end


function add_one_variable(poly, newparent)
    R = parent(poly)
    base = base_ring(R)

    t = last(gens(newparent))

    polybuilder = MPolyBuildCtx(newparent)
    for (e, c) in zip(exponent_vectors(poly), coefficients(poly))
        push_term!(polybuilder, c, [e..., 0])
    end

    return finish(polybuilder), t
end

function erase_last_variable(poly, newparent)
    # assuming poly is indepent of last variable

    R = parent(poly)
    base = base_ring(R)

    polybuilder = MPolyBuildCtx(newparent)
    for (e, c) in zip(exponent_vectors(poly), coefficients(poly))
        push_term!(polybuilder, c, e[1:end-1])
    end

    return finish(polybuilder)
end

# TODO : todo
function change_parent_ring(poly, newparent)
    original = parent(poly)
    originalnvars = nvars(original)
    
    parentnvars = nvars(newparent)
    
    if originalnvars == parentnvars
                
    elseif originalnvars + 1 == parentnvars
        return add_one_variable(poly, newparent) 
    elseif originalnvars == parentnvars + 1
        return erase_last_variable(poly, newparent)
    else
        @warn "failed to coerce polynomial $poly from $original to $parent"
    end
end


function unknown2known(u)
    libSingular.julia(libSingular.cast_number_to_void(u.ptr))    
end

function julia(x::Singular.n_FieldElem{Singular.FieldElemWrapper{Nemo.GaloisField,Nemo.gfp_elem}})
    libSingular.julia(libSingular.cast_number_to_void(x.ptr)).data
end

function julia(x::Singular.N_Field{Singular.FieldElemWrapper{U, T}}) where {U, T}
    libSingular.julia(x.ptr).data 
end

function julia(x)
    x
end


function singular2aa(poly::Singular.spoly{T}; base=false, new_ring=false) where {T}
    nvariables = length(gens(parent(poly)))
    xstrings = ["x$i" for i in 1:nvariables]
    if base == false
        base = base_ring(poly)
    end
    if new_ring == false
        new_ring, = AbstractAlgebra.PolynomialRing(base, xstrings, ordering=:lex)
    end
    change_base_ring(base, poly, parent=new_ring)
end

function double_singular2aa(
           poly::Singular.spoly{Singular.n_RingElem{Singular.spoly{T}}};
           base=false, new_ring=false) where {T}
    outer_change = singular2aa(poly)

    basebase = base_ring(parent(unknown2known(collect(coefficients(poly))[1])))

    nvariables = length(gens(parent(outer_change)))
    ystrings = ["y$i" for i in 1:nvariables]
    new_ring, = AbstractAlgebra.PolynomialRing(basebase, ystrings, ordering=:lex)

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
        new_ring, = Singular.PolynomialRing(base, xstrings, ordering=:lex)
    end
    change_base_ring(base, poly, parent=new_ring)    
end

function double_aa2singular(poly::MPoly{T}; base=false, new_ring=false) where {T}
    nvariables = length(gens(parent(poly)))
    ystrings = ["y$i" for i in 1:nvariables]
    xstrings = ["x$i" for i in 1:nvariables]
    basebase = base_ring(base_ring(parent(poly)))
    if base == false
        base, = Singular.PolynomialRing(basebase, xstrings, ordering=:lex)
    end
    if new_ring == false
        new_ring, = Singular.PolynomialRing(base, ystrings, ordering=:lex)
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

function applytofrac(f, x::Frac{T}) where {T}
    return map(f, (numerator(x), denominator(x)))
end

function tosingular(F::AbstractAlgebra.Generic.Rationals)
    Singular.QQ
end

function tosingular(F::Nemo.GaloisField)
    Singular.N_ZpField(Int(AbstractAlgebra.characteristic(F)))
end

function tosingular(F::Nemo.FlintRationalField)
    Singular.QQ
end

### TODO 
function tosingular(F::Nemo.GFPMPolyRing)
    strings = string.(gens(F))
    ground = base_ring(F)
    Singular.PolynomialRing(tosingular(ground), strings)[1]
end

function tosingular(F::AbstractAlgebra.Generic.MPolyRing)
    strings = string.(gens(F))
    ground = base_ring(F)
    Singular.PolynomialRing(tosingular(ground), strings)[1]
end

function tosingular(F)
    F
end

function toaa(::Singular.Rationals)
    AbstractAlgebra.Generic.QQ
end


###############################################################################

__randx = 10^4

function Base.rand(::AbstractAlgebra.Rationals{BigInt})
    return AbstractAlgebra.QQ(rand(1:__randx), rand(1:__randx))
end

function Base.rand(::Nemo.FlintRationalField)
    return Nemo.QQ(rand(1:__randx), rand(1:__randx))
end

function Base.rand(::Singular.Rationals)
    return Singular.QQ(rand(1:__randx), rand(1:__randx))
end


function Base.rand(wrapped::Singular.N_Field{Singular.FieldElemWrapper{U, T}}) where {U, T}
    return Base.rand( julia(wrapped) )
end

###############################################################################



