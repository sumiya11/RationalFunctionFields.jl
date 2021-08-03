
include("interpolation.jl")
include("utils.jl")
include("parse.jl")


import Singular: std, Ideal
import AbstractAlgebra
import Nemo
using DataStructures

import Singular: libSingular


function exponents_new_generating_set(genset)


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

    ybasering,  = AbstractAlgebra.PolynomialRing(ground, ystrings)
    n = 20

    # @info "" I

    ans = Dict()

    for (varidx, var) in enumerate(gens(basepolyring))

        gbs = []

        # generating interpolation points while fixing all variables but var
        points = [
            [ ground(i) for i in 1:nvariables ]
            for i in 1:n
        ]
        for i in 1:length(points)
            points[i][varidx] = ground(i)
        end
        
        lightring, = Singular.PolynomialRing(tosingular(ground), ystrings)

        for point in points
            # @info "point $point"

            Is = [
                  map_coefficients(c -> evaluate(c, point), f)
                  for f in I
            ]

            Is = [
                change_base_ring(tosingular(ground), f, parent=lightring)
                for f in Is
            ]

            # @info "" Is
            Ideal = Singular.Ideal(lightring, Is)

            # @info "" Ideal
            gb = Singular.std(Ideal, complete_reduction=true)

            push!(gbs, gb)
        end

        # @info "" gbs

        gbs = [
            [
                change_base_ring(ground, f, parent=basepolyring)
                for f in gens(gb)
            ]
            for gb in gbs
        ]

        # @info "" gbs

        ints = []

        # println( "gbs[1]\n", gbs[1] )

        for (j, g) in enumerate(gbs[1])
            push!(ints, Dict())

            lcoeff = leading_coefficient(g)
            g = map_coefficients(c -> c // lcoeff, g)

            exps = collect(exponent_vectors(g))
            for (k, c) in enumerate(monomials(g))
                # @info "interpolating $k th coeff in $j th gen"

                xs = [ point[varidx] for point in points ]
                ys = [ coeff(gb[j], k) for gb in gbs ]

                UNI, = AbstractAlgebra.PolynomialRing(ground, "x")
                f = interpolate_rational_function(UNI, xs, ys)

                # @info "interpolated " f
                ints[j][exps[k]] = f
            end
        end

        ans[var] = [ ints,  gbs[1] ]

    end

    return ans
end


function add_one_variable(polys)
    R = parent(polys[1])
    base = base_ring(R)

    strings = map(String, symbols(R))

    newring, vs = AbstractAlgebra.PolynomialRing(base, [strings..., "t"])
    t = last(vs)

    newpolys = []

    for f in polys
        polybuilder = MPolyBuildCtx(newring)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            push_term!(polybuilder, c, [e..., 0])
        end
        push!(newpolys, finish(polybuilder))
    end

    return newpolys, t
end

function erase_one_variable(newring, polys)
    # assuming polys are indepent of last variable

    R = parent(polys[1])
    base = base_ring(R)

    newpolys = []

    for f in polys
        polybuilder = MPolyBuildCtx(newring)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            push_term!(polybuilder, c, e[1:end-1])
        end
        push!(newpolys, finish(polybuilder))
    end

    return newpolys
end




function naive_new_generating_set(genset)

    #=
        for fi in polys
        Fi = Q*fi

        so that Fi is a polynomial
    =#

    originalring = parent(numerator( first(genset) ))
    nvariables = length(gens(originalring))
    xstrings = ["x$i" for i in 1:nvariables] 

    # reduce coeffs to GF(..)    
    
    basepolyring = originalring

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

    It, t = add_one_variable(I)
    sat = 1 - Q*t

    It = [It..., sat]

    basepolyrings, = Singular.AsEquivalentSingularPolynomialRing(basepolyring)
    yoverxs,  = Singular.PolynomialRing(basepolyrings, [ystrings..., "t"])

    Is = [
          map_coefficients(
              c -> change_base_ring(FF, c, parent=basepolyrings),
              f
          )
          for f in It
   ]

    Is = map(F -> change_base_ring(base_ring(yoverxs), F, parent=yoverxs), Is)

    Ideal = Singular.Ideal(yoverxs, Is...)

    gb = collect(gens(Singular.std(Ideal, complete_reduction=true)))

    gb = [
        double_singular2aa(f, base=basepolyring, new_ring=parent(t))
        for f in gb
    ]

    gb_no_t = []
    for f in gb
        if degree(f, t) == 0
            push!(gb_no_t, f)
        end
    end

    gb = erase_one_variable(yoverx, gb_no_t)
    
    generators = []
    for poly in gb
        for c in coefficients(poly)
            push!(generators, c)
        end
    end

    return gb, generators
end


function compare_degrees(true_basis, interpolated_basis)

    total = 0
    hit = 0
    
    
    for (varidx, var) in enumerate(keys(interpolated_basis))

        cffs, basis = interpolated_basis[var]

        true_coeffs = []


        for f in gens(true_basis)
            push!(true_coeffs, Dict())


            f = double_singular2aa(f)
            lcoeff = leading_coefficient(f)
            f = map_coefficients(c -> c // lcoeff, f)

            true_exps = [collect(exponent_vectors(m)) for m in monomials(f)]

            for (i, c) in enumerate(coeffs(f))
                true_coeffs[end][true_exps[i]...] = c
            end

        end

        Rxs = base_ring(true_coeffs[1][[0, 0]])
        xsvars = collect(gens(Rxs))
        
        ground = base_ring(Rxs)

        xstrings = ["x$i" for i in 1:length(xsvars)]
        Rxsaa, myxs = AbstractAlgebra.PolynomialRing(ground, xstrings)

        x = xsvars[varidx]
        myx = myxs[varidx]


        myeverything = Dict{Any, Set}()
        trueverything = Dict{Any, Set}()

        for (f, g) in zip(cffs, true_coeffs)
            for e in keys(f)
                if !haskey(myeverything, e)
                    myeverything[e] = Set()
                end
                if !haskey(trueverything, e)
                    trueverything[e] = Set()
                end

                myterm = f[e]

                trueterm = 1
                if !haskey(g, e)
                    for gg in true_coeffs
                        if haskey(gg, e)
                            trueterm = gg[e]
                            break
                        end
                    end
                else
                    trueterm = g[e]
                end

                mydeg = max(degree(numerator(myterm)),
                           degree(denominator(myterm)))


                xs = collect(gens(parent(numerator(trueterm))))

                truedeg = max(degree(numerator(trueterm), xs[varidx]),
                           degree(denominator(trueterm), xs[varidx]))

                push!(myeverything[e], mydeg)
                push!(trueverything[e], truedeg)
            end
        end


        for (key, val) in myeverything
            for dg in val
                if haskey(trueverything, key)

                    if dg in trueverything[key]
                        hit += 1
                    else
                        @warn "bedabedaogorchenie"
                        # return hit, total
                    end

                else
                    @warn "beda"
                    return hit, total
                end

                total += 1


            end
        end
    end


    return hit, total
end


function display_generating_exponents(genset)
    inter = exponents_new_generating_set(genset)

    answer = [ Dict() for _ in 1:length(first(values(inter))[1]) ]

    for v in keys(inter)
        exps, basis = inter[v]

        for (idx, g) in enumerate(exps)
            for ee in keys(g)
                if !haskey(answer[idx], ee)
                    answer[idx][ee] = 0
                end

                dg = degree(g[ee])

                answer[idx][ee] += dg

            end

        end

    end

    return answer
end



function new_generating_set(genset)

    #=
        TODO:
            . ensure evaluated generators are of the same structure
            . handle exponents better
            .
    =#

    basepolyring = parent(numerator( first(genset) ))
    nvariables = length(gens(basepolyring))
    
    ground = base_ring(basepolyring)

    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end

    Fs = map(numerator ∘ (g -> g * divexact(Q, denominator(g))), genset)

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

    I, t = add_one_variable(I)
    sat = 1 - Q*t

    I = [I..., sat]

    exponents = display_generating_exponents(genset)
    maxexp = maximum([maximum(values(s)) for s in exponents])

    gbs = []

    npoints = nvariables^(maxexp+1) + 2


    ndeg = maxexp
    nvs = length(gens(basepolyring))
    points = [ [ BigInt(j)^i for i in [ (ndeg+1)^k for k in 1:nvs ] ] for j in 1:npoints ]
    
    
    lightring, = Singular.PolynomialRing(tosingular(ground), [ystrings..., "t"])

    for point in points
        Is = [
                map_coefficients(c -> evaluate(c, point), f)
                for f in I
        ]

        Is = [
                change_base_ring(tosingular(ground), f, parent=lightring)
                for f in Is
        ]

        Ideal = Singular.Ideal(lightring, Is)

        gb = Singular.std(Ideal, complete_reduction=true)
        push!(gbs, collect(gens(gb)))
    end

    extended_ring, vvs= AbstractAlgebra.PolynomialRing(ground, [ystrings..., "t"])
    t = vvs[end]

    gbs = [
        [
           change_base_ring(ground, f, parent=extended_ring)
           for f in gb
        ]
        for gb in gbs
    ]

    gbs_no_t = []
    for gb in gbs
        gb_no_t = []
        for f in gb
            if degree(f, t) == 0
                push!(gb_no_t, f)
            end
        end
        push!(gbs_no_t, gb_no_t)
    end
    gbs = gbs_no_t

    gbs = [
            erase_one_variable(basepolyring, gb)
            for gb in gbs
          ]

    ints = []

    uni, = AbstractAlgebra.PolynomialRing(ground, "x")

    xs = []
    ys = []

    for i in 1:npoints
        p = points[i]

        push!(xs, p[1])
        push!(ys, [Dict() for _ in 1:length(gbs[i])])

        for (j, gg) in enumerate(gbs[i])
            for (ee, c) in zip(exponent_vectors(gg), coefficients(gg))
                ys[i][j][ee] = c
            end
        end

    end

    answer_e = [Dict() for _ in 1:length(gbs[1])]
    for (j, gg) in enumerate(gbs[1])

        for ev in exponent_vectors(gg)

            yssmall = [
                haskey(y[j], ev) ? y[j][ev] : ground(0)
                for y in ys
            ]


            f = interpolate_multivariate_rational_function(
                    basepolyring,
                    uni,
                    xs,
                    yssmall,
                    8)
            # TODO ??? or not to do: handle exponents better ? 
            

            R = basepolyring
            lol = [R(0), R(0)]
            for (i, D) in [ [1, numerator(f)], [2, denominator(f)] ]
                polybuilder = MPolyBuildCtx(R)
                for (e, c) in zip(0:degree(D), coefficients(D))
                    if iszero(c)
                        continue
                    end
                    push_term!(polybuilder, c, decompose_by_degrees(e, ndeg+1, nvs))
                end
                lol[i] = finish(polybuilder)
            end

            F = lol[1] // lol[2]
            answer_e[j][ev] = F

        end
    end

    gb = []
    R, = AbstractAlgebra.PolynomialRing(FractionField(basepolyring), ystrings)
    for etoc in answer_e
        polybuilder = MPolyBuildCtx(R)
        for e in keys(etoc)
            push_term!(polybuilder, etoc[e], e)
        end
        push!(gb, finish(polybuilder))
    end
    
    generators = []
    for poly in gb
        for c in coefficients(poly)
            push!(generators, c)
        end
    end

    return gb, generators
end















