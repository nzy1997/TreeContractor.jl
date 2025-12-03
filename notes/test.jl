using OMEinsum.OMEinsumContractionOrders, OMEinsum

code = ein"abc,bfg,ceh,cde->"
optcode = optimize_code(code, uniformsize(code, 2), OMEinsum.TreeSA())

label_elimination_order = OMEinsum.label_elimination_order(optcode)