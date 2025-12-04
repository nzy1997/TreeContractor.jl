using TreeContractor
using TreeContractor.OMEinsum
using ProblemReductions, Graphs

"""
    qubo_to_tensors(qubo::ProblemReductions.QUBO)

Convert a QUBO problem to tensor network representation.
For a QUBO: minimize ∑ᵢⱼ Qᵢⱼ xᵢ xⱼ
The QUBO matrix has diagonal elements representing linear terms and off-diagonal elements representing quadratic terms.
We create tensors where:
- Each variable xᵢ has a tensor with weight exp(-βE(xᵢ))
- Linear terms (diagonal) contribute to rank-1 tensors
- Quadratic terms (off-diagonal) contribute to rank-2 tensors
"""
function qubo_to_tensors(qubo::ProblemReductions.QUBO, β::Float64=1.0)
    Q = qubo.matrix
    n = size(Q, 1)
    
    # Create labels for each variable (1 to n)
    labels = collect(1:n)
    
    # Initialize tensors list
    tensors = Vector{Array{Float64}}()
    einsum_labels = Vector{Vector{Int}}()
    
    # Add linear terms (rank-1 tensors for each variable from diagonal)
    for i in 1:n
        # Create tensor for variable i with linear bias Q[i,i]
        # tensor[0] = exp(-β * Q[i,i] * 0) = 1
        # tensor[1] = exp(-β * Q[i,i] * 1) = exp(-β * Q[i,i])
        linear_tensor = [1.0, exp(-β * Q[i, i])]
        push!(tensors, linear_tensor)
        push!(einsum_labels, [i])
    end
    
    # Add quadratic terms (rank-2 tensors for off-diagonal elements)
    # Note: QUBO matrix is symmetric, so we only need upper triangle
    for i in 1:n
        for j in (i+1):n
            if Q[i, j] != 0.0  # Only add non-zero interactions
                # Create rank-2 tensor for interaction between variables i and j
                # tensor[xi, xj] = exp(-β * Q[i,j] * xi * xj)
                quad_tensor = zeros(2, 2)
                for xi in 0:1
                    for xj in 0:1
                        quad_tensor[xi+1, xj+1] = exp(-β * Q[i, j] * xi * xj)
                    end
                end
                
                push!(tensors, quad_tensor)
                push!(einsum_labels, [i, j])
            end
        end
    end
    
    return tensors, einsum_labels, labels
end

# Example: Create a fully connected QUBO problem
println("="^60)
println("Fully Connected QUBO Problem with Tensor Networks")
println("="^60)

# Create a complete graph with 20 nodes
n = 20
graph = complete_graph(n)

# Define QUBO coefficients
# Linear terms (h): bias for each variable
h = ones(n)

# Quadratic terms (Q): interaction weights for each edge
# For complete graph with n nodes, we have n*(n-1)/2 edges
Q = ones(n * (n - 1) ÷ 2)

qubo_problem = ProblemReductions.QUBO(graph, Q, h)

println("\nQUBO Problem:")
println("  Number of variables: $n")
println("  Number of linear terms: $(length(h))")
println("  Number of quadratic terms: $(length(Q))")
println("  Number of edges: $(ne(graph))")

# Step 1: Convert QUBO to tensor network
tensors, einsum_labels, labels = qubo_to_tensors(qubo_problem, 1.0)

# Step 2: Create einsum expression and optimize contraction order
ein_expr = EinCode(einsum_labels, Int[])
size_dict = OMEinsum.get_size_dict(einsum_labels, tensors)
optcode1 = optimize_code(ein_expr, size_dict, OMEinsum.TreeSA())
optcode2 = optimize_code(ein_expr, size_dict, OMEinsum.PathSA())

println("\nContraction complexity:")
println("TreeSA:\n", contraction_complexity(optcode1, size_dict))
println("PathSA:\n", contraction_complexity(optcode2, size_dict))

# Step 3: Contract the tensor network
Z1 = @time optcode1(tensors...)[]
Z2  = @time contract_with_mps(optcode2, tensors, size_dict; maxdim=100)[][]