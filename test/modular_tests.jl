

using .RationalFunctionFields: modular_reduction, rational_reconstruction


@testset "reduction + reconstruction tests" begin
    modulo = 2^31 - 1
    
    FF = Sing.GF(modulo)

    xs  = [ AA.QQ(rand(1:10000), rand(1:10000)) for _ in 1:10 ]
    xrs = modular_reduction(xs, FF)    
    ys  = rational_reconstruction(xrs, BigInt(modulo))
    
    @test xs == ys
    
    xs  = [ Sing.QQ(rand(1:10000), rand(1:10000)) for _ in 1:10 ]
    xrs = modular_reduction(xs, FF)
    ys  = rational_reconstruction(xrs, BigInt(modulo))

    @test xs == xrs    

    xs  = [ Nemo.QQ(rand(1:10000), rand(1:10000)) for _ in 1:10 ]
    xrs = modular_reduction(xs, FF) 
    ys  = rational_reconstruction(xrs, BigInt(modulo)) 
    
    @test xs == ys
end

@testset "complex reduction tests" begin
    R, (x, y) = AA.PolynomialRing(AA.QQ, ["x", "y"])
    
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




