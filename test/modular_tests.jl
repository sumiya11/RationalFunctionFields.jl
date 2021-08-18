

using .RationalFunctionFields: modular_reduction, rational_reconstruction, CRT,
                                CRT_reconstruction


@testset "reduction + reconstruction tests" begin
    modulo = 2^31 - 1
    
    FF = Sing.GF(modulo)

    xs  = [ Sing.QQ(rand(1:10000), rand(1:10000)) for _ in 1:10 ]
    xrs = modular_reduction(xs, FF)
    ys  = rational_reconstruction(xrs, BigInt(modulo))

    @test xs == ys    

end

@testset "complex reduction tests" begin
    R, (x, y) = AA.PolynomialRing(Sing.QQ, ["x", "y"])
    
    modulo = 2^31 - 1
    FF = AA.GF(2^31 - 1)

    fs = [
        R(0),
        R(1),
        7x^2*y + 5x + 3y,
        (x + y)^5
    ]
    
    frs = modular_reduction(fs, FF)
    
    fr = frs[3]
    @test AA.base_ring(fr) == FF
    
    frs = rational_reconstruction(frs, BigInt(modulo))
    
    @test fs == frs
end

@testset "CRT tests" begin

    moduli = [11]
    xs = [3, 4, 5, 10]

    for x in xs
        rems = [mod(x, mi) for mi in moduli]
        @test CRT(moduli, rems) == x
    end

    moduli = Primes.primes(5, 20)
    m = prod(moduli)
    
    xs = [ rand(1:(m-1)) for _ in 1:100 ]
    
    for x in xs
        rems = [mod(x, mi) for mi in moduli]
        
        @test CRT(moduli, rems) == x
    end

    moduli = Primes.primes(2^5, 2^6)
    m = prod(moduli)

    xs = [ rand(1:(m-1)) for _ in 1:100 ]

    for x in xs
        rems = [mod(x, mi) for mi in moduli]

        @test CRT(moduli, rems) == x
    end

end


@testset "CRT basis reconstruction" begin
    moduli = [7, 5]
    
    FF1 = Sing.N_ZpField(moduli[1])
    FF2 = Sing.N_ZpField(moduli[2])
    
    R1, (x, y) = AA.PolynomialRing(FF1, ["x1", "y1"])
    R11, (a, b) = AA.PolynomialRing(AA.FractionField(R1), ["a1", "b1"])
    
    fs1 = [
        4*a^2*b^2 + (x + 6y)//(x + y)*a + (3x + 2x*y + 1*y)//(y)*b + 5y^2,
        1a + 3b + 2a*b + 2x//(1y)
    ]
    
    R2, (x, y) = AA.PolynomialRing(FF2, ["x2", "y2"])
    R22, (a, b) = AA.PolynomialRing(AA.FractionField(R2), ["a2", "b2"])

    fs2 = [
        2*a^2*b^2 + (1x + 3y)//(x + y)*a + (3x + 1x*y + 4*y)//(y)*b + 5y^2,
        1a + 3b + 1a*b + 4x//(y) 
    ]

    basis = CRT_reconstruction([fs1, fs2], moduli) 
    
    yoverx = AA.parent(basis[1])
    basepolyrings = AA.base_ring(yoverx)
    ground = AA.ZZ

    x, y = AA.gens(AA.base_ring(basepolyrings))
    a, b = AA.gens(yoverx)

    # to check when number is in the denominator
    ans = [
           32*a^2*b^2 + (x + 13*y)//(x + y)*a + (16*x*y + 3*x + 29*y)//(y)*b + 5*y^2,
           1a + 3b + 16a*b + 9x//(y)
    ]

    @test basis == ans

end































