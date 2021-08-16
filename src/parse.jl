
"""
    Loads a set of generators from the file by given filepath
    One could check valid file examples in the RFF/data directory
    
    Returns a set of fractions of AA polys over AA rationals
"""
function load_generators(filepath)
    lines = []
    open(filepath, "r") do inputs
        lines = map(strip, readlines(inputs))
    end

    @info "Loading $filepath"
    
    strings = map(String, split(lines[1], ", "))

    S, xs = AbstractAlgebra.PolynomialRing(Singular.QQ, strings, ordering=:lex)
    
    mapping = Dict{Symbol, AbstractAlgebra.Generic.MPoly}(
        Symbol(x) => AbstractAlgebra.gen(S, i)
        for (i, x) in enumerate(strings)
    )
    
    generators = []

    for line in lines[2:end]
        polystrings = map(Meta.parse ∘ String ∘ strip, split(line, ", "))
        polys = map(f -> myeval(f, mapping), polystrings)
        
        append!(generators, [ g // S(polys[1]) for g in polys[2:end] ])
    end
       
    @info "loaded $(length(generators)) from $filepath"
    
    return generators
end





















