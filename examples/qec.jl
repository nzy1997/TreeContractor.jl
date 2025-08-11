using TreeContractor
using TensorQEC
using TreeContractor.OMEinsum
using Random

dem = TensorQEC.parse_dem_file("examples/data/surface_code_d=3_r=3.dem")

ct = compile(TNMMAP(TreeSA(),true),dem)

contraction_complexity(ct.optcode,uniformsize(ct.optcode,2))

# Random.seed!(1234)
# ep = random_error_qubits(dem)
# syd = syndrome_extraction(ep, ct.tanner)

# TensorQEC.update_syndrome!(ct,syd)
# TensorQEC.update_syndrome!(ct2,syd)

@time ct.optcode(ct.tensors...)
#   0.007278 seconds (113.39 k allocations: 7.312 MiB)
# 2-element Vector{Float64}:
#  0.8426216150912986
#  9.535387828248582e-8

ct2 = compile(TNMMAP(OMEinsum.PathSA(),true),dem)

@time contract_with_mps(ct2.optcode, ct2.tensors, uniformsize(ct2.optcode,2); maxdim = 20)
@time contract_with_mps(ct2.optcode, ct2.tensors, uniformsize(ct2.optcode,2); maxdim = 50)
@time contract_with_mps(ct2.optcode, ct2.tensors, uniformsize(ct2.optcode,2); maxdim = 80)