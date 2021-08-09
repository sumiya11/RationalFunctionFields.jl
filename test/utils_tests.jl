

if !isdefined(Main, :testset)
    using Test
    using TestSetExtensions
end

include("../src/RationalFunctionFields.jl")
using .RationalFunctionFields: change_parent_ring, singular2aa, aa2singular, 
                        double_singular2aa, double_aa2singular


import Singular
import AbstractAlgebra

@testset "Utils tests" begin

    R, xs = Singular.PolynomialRing(Singular.QQ, ["x", "y"])
    x, y = xs
    f = x + y

    faa = singular2aa(f)

    revfaa = aa2singular(faa, new_ring=R)
    @test revfaa == f    

    R2, xs2 = Singular.PolynomialRing(R, ["x2", "y2"])
    x2, y2 = xs2
    f2 = x2 + y2

    f2aa = double_singular2aa(f2)

    revf2s = double_aa2singular(f2aa, base=R, new_ring=R2)
    @test revf2s == f2

    
    AA, (x, y) = AbstractAlgebra.PolynomialRing(Singular.QQ, ["x", "y"])
    f = x + y
    newparent, (xx, yy, tt) = AbstractAlgebra.PolynomialRing(Singular.QQ, ["x", "y", "t"])
    
    ft, t = change_parent_ring(f, newparent)
    
    @test parent(ft) == newparent 
    @test t == tt

    frev = change_parent_ring(ft, AA)

    @test f == frev


end



