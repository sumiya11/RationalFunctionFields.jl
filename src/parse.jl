
"""
    Loads a set of generators from the file by given filepath
    One could check valid file examples in the RFF/data directory
    
    Returns a set of fractions of Nemo polys over Nemo rationals
"""
function load_generators(filepath)
    lines = []
    open(filepath, "r") do inputs
        lines = map(strip, readlines(inputs))
    end

    @info "Loading from $filepath"
    
    strings = map(String, split(lines[1], ", "))

    S, xs = AbstractAlgebra.PolynomialRing(AbstractAlgebra.QQ, strings, ordering=:lex)

    mapping = Dict{Symbol, MPoly}(
        Symbol(x) => AbstractAlgebra.gen(S, i)
        for (i, x) in enumerate(strings)
    )
    
    nemoring, = Nemo.PolynomialRing(Nemo.QQ, strings, ordering=:lex)
    
    generators = []

    for line in lines[2:end]
        polystrings = map(Meta.parse ∘ String ∘ strip, split(line, ", "))
        polys = map(f -> myeval(f, mapping), polystrings)
       
        polys = [
            Nemo.isconstant(S(f)) ? nemoring(Nemo.QQ(f)) : change_base_ring(Nemo.QQ, f, parent=nemoring)
            for f in polys
        ]

        append!(generators, [ g // polys[1] for g in polys[2:end] ])
    end
       
    @info "loaded $(length(generators)) generators from $filepath"
    
    return generators
end





















