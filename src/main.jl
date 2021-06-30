
include("interpolation.jl")


using Singular
using AbstractAlgebra
using Nemo


GROUND = Singular.QQ


function exponents_new_generating_set(genset)

    basepolyring = parent(numerator( first(genset) ))
    nvariables = length(gens(basepolyring))

    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end

    Fs = map(numerator ∘ (g -> g * Q), genset)

    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = AbstractAlgebra.PolynomialRing(basepolyring, ystrings)

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    I = [
        Fyi*Qx - Fxi*Qy
        for (Fyi, Fxi) in zip(Fx, Fy)
    ]
    
    # substitute + groebner n times
    
    ybasering,  = AbstractAlgebra.PolynomialRing(GROUND, ystrings)
    n = 20
    
    @info "" I
        
    ans = Dict()
    
    for (varidx, var) in enumerate(gens(basepolyring))
        @info "$varidx th, handing $var"
        
        gbs = []
        
        # generating interpolation points while fixing all variables but var
        points = [
            [ GROUND(0) for _ in  1:nvariables ]
            for _ in 1:n
        ]
        for i in 1:length(points)
            points[i][varidx] = GROUND(i)
        end
        
        lightring, = Singular.PolynomialRing(GROUND, ystrings)
        
        for point in points
            @info "point $point"
    
            Is = [
                  map_coefficients(c -> evaluate(c, point), f)
                  for f in I
            ]
            
            Is = [
                change_base_ring(GROUND, f, parent=lightring) 
                for f in Is
            ]
            
            @info "" Is
            Ideal = Singular.Ideal(lightring, Is)

            @info "" Ideal
            gb = std(Ideal)

            push!(gbs, gb)
        end
        
        @info "" gbs
        
        gbs = [
            [
                change_base_ring(GROUND, f, parent=basepolyring)
                for f in gens(gb)
            ]
            for gb in gbs
        ]

        @info "" gbs
        
        ints = []
        
        println( "gbs[1]\n", gbs[1] )

        for (j, g) in enumerate(gbs[1])
            push!(ints, [])
            for (k, c) in enumerate(coeffs(g))
                @info "interpolating $k th coeff in $j th gen"
                
                xs = [ point[varidx] for point in points ]
                ys = [ coeff(gb[j], k) for gb in gbs ]
                
                UNI, = AbstractAlgebra.PolynomialRing(GROUND, "x")
                f = interpolate_rational_function(UNI, xs, ys)
                
                @info "interpolated " f
                
                push!(ints[j], f)
            end
        end
        
        ans[var] = [ ints, [ gbs[1], gbs[2], gbs[3] ] ]
        
    end
    
    return ans
end



function naive_new_generating_set(genset)
    
    #=
        for fi in polys
        Fi = Q*fi

        so that Fi is a polynomial
    =#

    basepolyring = parent(numerator( first(genset) ))
    nvariables = length(gens(basepolyring))
    
    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end

    Fs = map(numerator ∘ (g -> g * Q), genset)

    #=
        I = <
            Fi(y) Q(x) - Fi(x) Q(y)
        >
    =#

    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = AbstractAlgebra.PolynomialRing(basepolyring, ystrings)

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    I = [
        Fyi*Qx - Fxi*Qy
        for (Fyi, Fxi) in zip(Fx, Fy)
    ]

    basepolyrings, = Singular.AsEquivalentSingularPolynomialRing(basepolyring)
    yoverxs,  = Singular.PolynomialRing(basepolyrings, ystrings) 
    
    Is = [
          map_coefficients(
                    c -> change_base_ring(GROUND, c, parent=basepolyrings),
                    f
          )
          for f in I  
   ]
    
    Is = map(F -> change_base_ring(base_ring(yoverxs), F, parent=yoverxs), Is)

    Ideal = Singular.Ideal(yoverxs, Is...)

    gb = std(Ideal)

    return gb
end

R, (x1, x2) = AbstractAlgebra.PolynomialRing(GROUND, ["x1", "x2"])

Q = FractionField(R)

f1 = (x1 + x2^2) // 1
f2 = (x1 - x2) // 1

exps = exponents_new_generating_set([f1, f2])
ggen = naive_new_generating_set([f1, f2])


for (var, smth) in exps
    println(var)
    println("coeffs = ", smth[1])
    println("possible bases = \n", join(map(string, smth[2]), '\n'))
end

println( "\ntrue\n", ggen )
























