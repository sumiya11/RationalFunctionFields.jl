


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
    evalring    

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
    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = Nemo.PolynomialRing(
                                         basepolyring,
                                         ystrings,
                                         ordering=:lex
    )

    context = IdealContext(
                         [],
                         yoverx,
                         basepolyring,
                         nvariables,
                         ground,
                         ystrings,
                         1,
                         saturated,
                         1,
                         1
    )

    idealize(genset, context, saturated=saturated)
end


function idealize(
               genset::Array{T}, context::IdealContext;
               saturated=false) where {T}

    yoverx = context.yoverx
    basepolyring = context.basepolyring
    nvariables = context.nvariables
    ground = context.ground
    ystrings = context.ystrings

    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end

    Fs = map(numerator ∘ (g -> g * Q), genset)

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    I = [
        Fyi * Qx - Fxi * Qy
        for (Fyi, Fxi) in zip(Fy, Fx)
    ]
    
    # TODO: align saturation types    
    if saturated
        I, t = saturate(I, Q)
    else
        t = Nemo.gen(yoverx, 1)
    end
    
    evalstrings = deepcopy(ystrings)
    if saturated
        push!(evalstrings, "t")
    end

    singular_ground = tosingular(ground)
    evalring, evalvars = Singular.PolynomialRing(
                       singular_ground,
                       evalstrings,
                       ordering=ordering_lp(nvariables)*ordering_c()
    )

    return IdealContext(
               I,
               yoverx,
               basepolyring,
               nvariables,
               ground,
               ystrings,
               evalring,
               saturated,
               t,
               Q
    )
end

function idealize(
                elem, context::IdealContext;
                saturated=false) where {T}
    idealize([elem], context)
end























