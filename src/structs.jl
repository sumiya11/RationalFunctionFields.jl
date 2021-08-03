
include("main.jl")


###############################################################################


mutable struct RationalFunctionField
    generating_set     # present always
    groebner_coeffs    # lazy field, not to be addressed directly
    groebner_ideal     # lazy field, not to be addressed directly
end


function RationalFunctionField(genset)
    RationalFunctionField(genset, [], [])
end

function gens(FF::RationalFunctionField)
    return FF.groebner_coeffs
end

function contains_randomized(FF::RationalFunctionField, elem)
    # pass
end


function contains_using_groebner(FF::RationalFunctionField, elem)
    # elem = A // B
    A, B = numerator(elem), denominator(elem)
    
    I = FF.groebner_ideal
    
    nvs = nvars(parent(A))
    xstrings = ["x$i" for i in 1:nvs]
    ystrings = ["y$i" for i in 1:nvs]
   
    basepolyring = parent(A)
    F = base_ring(basepolyring)
    
    yoverx, yoverxvars = AbstractAlgebra.PolynomialRing(basepolyring, ystrings)

    Ay = change_base_ring(basepolyring, A, parent=yoverx)
    By = change_base_ring(basepolyring, B, parent=yoverx)

    basepolyrings, = Singular.AsEquivalentSingularPolynomialRing(basepolyring)
    fracbases = FractionField(basepolyrings)
    yoverxs, = Singular.PolynomialRing(fracbases, ystrings)
    
    Is = [
          map_coefficients(c -> change_base_ring(F, numerator(c), parent=basepolyrings) // change_base_ring(F, denominator(c), parent=basepolyrings), f) for f in I 
    ]

    Is = map(c -> change_base_ring(base_ring(yoverxs), c, parent=yoverxs), Is)

    f = change_base_ring(base_ring(yoverxs), A * By - Ay * B, parent=yoverxs)  
    
    fs = Ideal(yoverxs, f)
    Is = Ideal(yoverxs, Is)
    
    return Singular.contains(Is, fs)
end


function compute_groebner!(
                               FF::RationalFunctionField;
                               backend_algorithm=new_generating_set)
    
    FF.groebner_ideal, FF.groebner_coeffs = backend_algorithm(FF.generating_set)
end


function simple_preprocess(FF::RationalFunctionField)
    generators = FF.groebner_coeffs
    tmp = copy(generators)
    tmp = filter((c) -> !isconstant(c), tmp)
    sort!(tmp, by=(c) -> sum(degrees(c)), rev=true)
    ans = [ first(tmp) ]
    for elem in tmp
        field = RationalFunctionField(ans)
        compute_groebner!(field)
        if !contains(field, elem)
            push!(ans, elem)
        end
    end
    ans
end

function simplify_generators!(
                              FF::RationalFunctionField;
                              backend_algorithm=simple_preprocess)
    FF.groebner_coeffs = backend_algorithm(FF)
end


function contains(FF::RationalFunctionField, elem; proved=true)
    if proved
        contains_using_groebner(FF, elem)
    else
        contains_randomized(FF, elem)
    end
end


###############################################################################

















