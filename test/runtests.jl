using TreeContractor
using Test
using OMEinsum
using Random
Random.seed!(1234)

code = ein"abc,cde,egh,bfg->"
size_dict = uniformsize(code, 2)
optcode = optimize_code(code, size_dict, TreeSA())


t1 = rand(2,2,2)
t2 = rand(2,2,2)
t3 = rand(2,2,2)
t4 = rand(2,2,2)

optcode(t1, t2, t3, t4)
code(t1, t2, t3, t4)

TreeContractor.generate_tree(optcode)