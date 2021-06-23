

using AbstractAlgebra

#=
function substitute_coeffs(mpoly, point, newbasering)

    # reducio & reconstruction
    for trms in terms(mpoly)
        # this is horrible, change
        c, e = coeff(trms, 1), exponent(trms, 1)
        cevaluated = evaluate(c, point)
        setcoeff!(mpoly, e , cevaluated)
    end

    return change_base_ring(QQ, mpoly, parent=newbasering)
end
=#


function new_generating_set(genset)

    #=
        for fi in polys
        Fi = Q*fi

        so that Fi is a polynomial
    =#

    basepolyring = parent(numerator( first(genset) ))
    nvariables = length(gens(basepolyring))

    Q = lcm(map(denominator, genset)...)
    Fs = map(numerator ∘ (g -> g * Q), genset)

    @info "" Q
    @info "" Fs

    #=
        I = <
            Fi(y) Q(x) - Fi(x) Q(y)
        >
    =#

    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = PolynomialRing(basepolyring, ystrings)

    Fx = map(yoverx, Fs)
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = yoverx(Q)
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    @info "" Fx
    @info "" Fy


    @info "" Qx
    @info "" Qy

    I = [
        Fyi*Qx - Fxi*Qy
        for (Fyi, Fxi) in zip(Fx, Fy)
    ]


    @info "" I

    #=
    lightpoly, yvars = PolynomialRing(QQ, ystrings)

    for evalpoint in evalpoints
    evaluatedI = [
        substitute_coeffs(poly, point, lightpoly)
        for poly in I
    ]
    =#

    # groebnerbasis = basis(evaluatedI)

    # interpolate

end

R, (x1, x2) = PolynomialRing(QQ, ["x1", "x2"])

Q = FractionField(R)

f1 = (x1 + x2) // x1
f2 = 1 // x2^2

new_generating_set([f1, f2])
