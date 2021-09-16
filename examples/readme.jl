
# The file contains example from readme

# -----------------------------------------------------------------------------

import Pkg;
Pkg.add(url="https://github.com/sumiya11/RationalFunctionFields.jl")

# -----------------------------------------------------------------------------

using Nemo
using RationalFunctionFields

# -----------------------------------------------------------------------------

R, (x, y) = PolynomialRing(QQ, ["x", "y"])

generators = [
            (x^2 + y^2) // R(1),
            1 // (x^2)
]

# -----------------------------------------------------------------------------

rff = RationalFunctionField(generators)

print( contains(rff, y^2, p=0.999) )

print( contains(rff, y^2, p=1) )














