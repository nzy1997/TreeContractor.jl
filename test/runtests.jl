using TreeContractor
using Test
using TreeContractor.OMEinsum
using Random
using LinearAlgebra


@testset "tensor2mps and contract_mps" begin
    ashape = (6, 3, 4, 5)
    a = rand(Float64, ashape)
    b, bd_vec = TreeContractor.tensor2mps(a)
    bd = 1
    @test bd_vec[1] == 1
    for (i,s) in enumerate(ashape)
        bdnew = size(b[i])[3]
        @test size(b[i]) == (bd, s, bdnew)
        bd = bdnew
        @test bd_vec[i+1] == bdnew
    end
    @test TreeContractor.contract_mps(b)[] ≈ sum(a) atol = 1e-10
end

@testset "delta_mps" begin
    a = TreeContractor.delta_mps(2,3,Float64)
    @test a[:,1,:] == I(2)
    @test a[:,2,:] == I(2)
    @test a[:,3,:] == I(2)
end

@testset "apply_rank_3_tensor" begin
    Random.seed!(1234)
    t1 = rand(2,2,2)
    t2 = rand(2,2,2)
    t3 = rand(2,2,2)
    t4 = rand(2,2,2)

    m = ein"bfc,efg,cad,gah->bedh"(t1,t2,t3,t4)
    m = reshape(m, 2*2,2*2)

    m1 = TreeContractor.apply_rank_3_tensor(t1,t2)
    m2 = TreeContractor.apply_rank_3_tensor(t3,t4)
    mp = ein"bfc,cad-> bd"(m1,m2)
    @test m ≈ mp atol = 1e-10
end

@testset "contract_with_mps" begin
    code = ein"abc,abd->"
    optcode = optimize_code(code, uniformsize(code, 2), OMEinsum.PathSA())

    Random.seed!(1234)
    t1 = rand(2,2,2)
    t2 = rand(2,2,2)

    tensors = [t1, t2]
    right_answer = optcode(t1, t2)[]

    @test contract_with_mps(optcode, tensors, uniformsize(code, 2))[1][] ≈ right_answer atol = 1e-10
end

@testset "contract_with_mps" begin
    code = ein"ac,bd,ab,cd->"
    size_dict = Dict('a' => 2, 'b' => 3, 'c' => 4, 'd' => 5)
    optcode = optimize_code(code, size_dict, OMEinsum.PathSA())

    Random.seed!(1234)
    t1 = rand(2,4)
    t2 = rand(3,5)
    t3 = rand(2,3)
    t4 = rand(4,5)

    tensors = [t1, t2, t3, t4]
    right_answer = optcode(tensors...)[]

    @test contract_with_mps(optcode, tensors, size_dict)[1][] ≈ right_answer atol = 1e-10
end

@testset "contract_with_mps" begin
    code = ein"abc,cde,egh,fbg->"
    optcode = optimize_code(code, uniformsize(code, 2), OMEinsum.PathSA())

    Random.seed!(1234)
    t1 = rand(2,2,2)
    t2 = rand(2,2,2)
    t3 = rand(2,2,2)
    t4 = rand(2,2,2)
    
    tensors = [t1, t2, t3, t4]
    right_answer = optcode(tensors...)[]
    @show right_answer

    mps, apply_vec, tensor_labels, vanish_labels_vec = TreeContractor.code2mps(optcode,uniformsize(code, 2)); mps = TreeContractor.apply_tensors!(mps, apply_vec, tensors, tensor_labels, vanish_labels_vec)

    @show mps.tensors

    @test contract_with_mps(optcode, tensors, uniformsize(code, 2))[1][] ≈ right_answer atol = 1e-10
end


@testset "canonicalize" begin
    Random.seed!(42)
    N, center, χ = 10, 5, 20
    T = ComplexF64
    mps = TreeContractor.random_mps(T, N; maxdim=χ)
    norm_init = norm(mps)

    TreeContractor.canonicalize!(mps, center)

    for i in 1:(center - 1)
        U = reshape(mps.tensors[i], size(mps.tensors[i])[1] * size(mps.tensors[i])[2], :)
        @test U' * U ≈ LinearAlgebra.I
    end
    for i in (center + 1):TreeContractor.nsite(mps)
        V = reshape(mps.tensors[i], :, size(mps.tensors[i])[2] * size(mps.tensors[i])[3])
        @test V * V' ≈ LinearAlgebra.I
    end

    @test TreeContractor.orthocenter(mps) == center
    @test norm(mps) ≈ norm_init
    @test TreeContractor.is_canonicalized(mps)
    @test TreeContractor.check_canonical(mps)
end


@testset "compression by local SVD" begin
    N = 6
    Random.seed!(42)
    χ = 8
    # Create MPS with high bond dimensions
    mps = TreeContractor.LabeledMPS([
        randn(ComplexF64, 1, 2, χ),
        [randn(ComplexF64, χ, 2, χ) for _ in 2:(N-1)]...,
        randn(ComplexF64, χ, 2, 1)
    ], [i for i in 1:N])
    
    original_vec = vec(mps)
    original_elements = TreeContractor.num_of_elements(mps)
    
    # Test compression
    TreeContractor.compress!(TreeContractor.LocalCompress(), mps; niters=2, maxdim=8)
    compressed_elements = TreeContractor.num_of_elements(mps)
    
    @test compressed_elements < original_elements
    @test TreeContractor.check_canonical(mps)
    @test vec(mps) ≈ original_vec atol=1e-4
end

@testset "compression by global SVD" begin
    N = 6
    Random.seed!(42)
    χ = 8
    # Create MPS with high bond dimensions
    mps = TreeContractor.LabeledMPS([
        randn(ComplexF64, 1, 2, χ),
        [randn(ComplexF64, χ, 2, χ) for _ in 2:(N-1)]...,
        randn(ComplexF64, χ, 2, 1)
    ], [i for i in 1:N])

    TreeContractor.normalize!(mps)
    
    original_vec = vec(mps)
    original_elements = TreeContractor.num_of_elements(mps)
    
    # Test compression
    TreeContractor.compress!(TreeContractor.FullCompress(), mps; maxdim=8)
    compressed_elements = TreeContractor.num_of_elements(mps)
    
    @test compressed_elements < original_elements
    @test TreeContractor.check_canonical(mps)
    @test vec(mps) ≈ original_vec atol=1e-14
end

@testset "contract_with_mps" begin
    code = ein"abc,cde,egh,fbg,yczd,ybzc->"
    bd = 3
    optcode = optimize_code(code, uniformsize(code, bd), OMEinsum.PathSA())

    Random.seed!(1234)
    t1 = rand(bd,bd,bd)
    t2 = rand(bd,bd,bd)
    t3 = rand(bd,bd,bd)
    t4 = rand(bd,bd,bd)
    t5 = rand(bd,bd,bd,bd)
    t6 = rand(bd,bd,bd,bd)

    tensors = [t1, t2, t3, t4, t5, t6]
    right_answer = optcode(tensors...)[]
    @show right_answer

    @test contract_with_mps(optcode, tensors, uniformsize(code, bd); maxdim = 10)[1][] ≈ right_answer atol = 1e-10
end


@testset "tensor with output label" begin
    code = ein"abc,cde,egh,fbg->f"
    bd = 5
    optcode = optimize_code(code, uniformsize(code, bd), OMEinsum.PathSA())

    Random.seed!(1234)
    t1 = rand(bd,bd,bd)
    t2 = rand(bd,bd,bd)
    t3 = rand(bd,bd,bd)
    t4 = rand(bd,bd,bd)
    
    tensors = [t1, t2, t3, t4]
    right_answer = optcode(tensors...)
    @show right_answer

    @test contract_with_mps(optcode, tensors, uniformsize(code, bd); maxdim = 10)[1][1,:,1] ≈ right_answer atol = 1e-10
end