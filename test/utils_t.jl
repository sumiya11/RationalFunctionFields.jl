
include("../src/utils.jl")


R, xs = Singular.PolynomialRing(Singular.QQ, ["x", "y"])
x, y = xs
f = x + y

faa = singular2aa(f)

revfaa = aa2singular(faa, new_ring=R)
@assert revfaa == f

R2, xs2 = Singular.PolynomialRing(R, ["x2", "y2"])
x2, y2 = xs2
f2 = x2 + y2

f2aa = deux_singular2aa(f2)

revf2s = deux_aa2singular(f2aa, base=R, new_ring=R2)
@assert revf2s == f2

