using TreeContractor
using Test
using OMEinsum
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

    @test contract_with_mps(optcode, tensors, uniformsize(code, 2)) ≈ right_answer atol = 1e-10
end

@testset "contract_with_mps" begin
    code = ein"ac,bd,ab->"
    size_dict = Dict('a' => 2, 'b' => 3, 'c' => 4, 'd' => 5)
    optcode = optimize_code(code, size_dict, OMEinsum.PathSA())

    Random.seed!(1234)
    t1 = rand(2,4)
    t2 = rand(3,5)
    t3 = rand(2,3)

    tensors = [t1, t2, t3]
    right_answer = optcode(tensors...)[]

    @test contract_with_mps(optcode, tensors, size_dict) ≈ right_answer atol = 1e-10
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
    # 11.468356722794542
    @show right_answer

    mps, apply_vec, tensor_labels = TreeContractor.code2mps(optcode,uniformsize(code, 2)); TreeContractor.apply_tensors!(mps, apply_vec, tensors, tensor_labels)

    @show TreeContractor.contract_mps(mps.tensors)

    @test contract_with_mps(optcode, tensors, uniformsize(code, 2)) ≈ right_answer atol = 1e-10
end

