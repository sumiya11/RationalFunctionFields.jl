# Adapted from https://discourse.julialang.org/t/expression-parser/41880/7
# code by Alan R. Rogers, Professor of Anthropology, University of Utah

# Adapted from https://github.com/x3042/Exact-reduction-of-ODE-systems/blob/main/src/myeval.jl
# by Alexander D. , Higher School of Economics

function myeval(e::Union{Expr,Symbol,Number}, map::Dict{Symbol,MPoly})
    try
        return _myeval(e, map)
    catch ex
        println("Can't parse \"$e\"")
        rethrow(ex)
    end
end

function _myeval(s::Symbol, map::Dict{Symbol,MPoly})
    if haskey(map, s)
        return map[s]
    else
        throw(UndefVarError(s))
    end
end

# Numbers are converted to QQ.
function _myeval(x::Number, map::Dict{Symbol,MPoly})
    return QQ(x)
end

# a helper definition for floats
function _myeval(x::Float64, map::Dict{Symbol, MPoly})
    result = QQ(0)

    # Getting the result from the string representation in order
    # to avoid approximations caused by the float representation
    s = string(x)
    denom = 1
    extra_num = 1
    if occursin(r"[eE]", s)
        s, exp = split(s, r"[eE]")
        if exp[1] == "-"
            denom = QQ(10)^(-parse(Int, exp))
        else
            extra_num = QQ(10)^(parse(Int, exp))
        end
    end
    frac = split(s, ".")
    if length(frac) == 1
        result = QQ(parse(fmpz, s)) * extra_num // denom
    else 
        result = QQ(parse(fmpz, frac[1] * frac[2])) * extra_num // (denom * 10^(length(frac[2])))
    end
    
    # too verbose for now
    # @warn "a possibility of inexact float conversion" from=x to=result
    return result
end

# To parse an expression, convert the head to a singleton
# type, so that Julia can dispatch on that type.
function _myeval(e::Expr, map::Dict{Symbol,MPoly})
    return _myeval(Val(e.head), e.args, map)
end

# Call the function named in args[1]
function _myeval(::Val{:call}, args, map::Dict{Symbol,MPoly})
    return _myeval(Val(args[1]), args[2:end], map)
end

# Addition
function _myeval(::Val{:+}, args, map::Dict{Symbol,MPoly})
    x = 0
    for arg ∈ args
        x += _myeval(arg, map)
    end
    return x
end

# Subtraction and negation
function _myeval(::Val{:-}, args, map::Dict{Symbol,MPoly})
    len = length(args)
    if len == 1
        return -_myeval(args[1], map)
    else
        return _myeval(args[1], map) - _myeval(args[2], map)
    end
end

# Multiplication
function _myeval(::Val{:*}, args, map::Dict{Symbol,MPoly})
    x = 1
    for arg ∈ args
        x *= _myeval(arg, map)
    end
    return x
end

# Division
function _myeval(::Val{:/}, args, map::Dict{Symbol,MPoly})
    return _myeval(args[1], map) // _myeval(args[2], map)
end

# Exponentiation
function _myeval(::Val{:^}, args, map::Dict{Symbol, MPoly})
    if 1 != denominator(_myeval(args[2], map))
        @warn "chto-to strannoe happened"
    end

    return _myeval(args[1], map) ^ Int(numerator(_myeval(args[2], map)))
end



