
using AbstractAlgebra


# univariate
# returns :: AbstractAlgebra.Generic.Poly
function interpolate(R, xs, ys)

end

# multivariate
# returns :: AbstractAlgebra.Generic.MPoly
function interpolate(R, xs, ys)

end


U, X = QQ["x"]
g = X^2

R, (x, y) = QQ["x", "y"]

f = x^2 + y
