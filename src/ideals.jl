


"""
    The structure represents the good ideal

"""
mutable struct IdealContext
    I
    yoverx
    basepolyring
    nvariables
    ground
    ystrings
    
    # if the the ideal is saturated
    saturated::Bool
    t
    Q
end


function idealize(elem)
    idealize([elem])
end

function idealize(genset::Array{T}; saturated=false) where {T}
    basepolyring = parent(numerator( first(genset) ))
    nvariables = length(gens(basepolyring))
    ground = base_ring(basepolyring)

    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end

    Fs = map(numerator ∘ (g -> g * Q), genset)

    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = Nemo.PolynomialRing(
                                         basepolyring,
                                         ystrings,
                                         ordering=:lex
    )

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    # Gleb: to cancel Fyi and Qy by their gcdi
    I = [
        Fyi * Qx - Fxi * Qy
        for (Fyi, Fxi) in zip(Fy, Fx)
    ]
    
    if saturated
        I, t = saturate(I, Q)
    else
        t = Nemo.gen(yoverx, 1)
    end

    return IdealContext(
               I,
               yoverx,
               basepolyring,
               nvariables,
               ground,
               ystrings,
               saturated,
               t,
               Q
    )
end






















