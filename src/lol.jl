
using Singular
using AbstractAlgebra

include("utils.jl")

GROUND = Singular.QQ

R, (x1, x2) = AbstractAlgebra.PolynomialRing(GROUND, ["x1", "x2"])

Q = FractionField(R)

f1 = (x1 + x2) // x1
f2 = 1 // x2^2

new_generating_set([f1, f2])










