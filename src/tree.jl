struct TensorBinaryTree{T,LT}
    tensor::Array{T, 3}
    left::Union{TensorBinaryTree, Nothing}
    right::Union{TensorBinaryTree, Nothing}
    label::LT
end

Base.show(io::IO, ::MIME"text/plain", p::TensorBinaryTree) = show(io, p)
function Base.show(io::IO, p::TensorBinaryTree)
    println(io, "label: $(p.label)")
    if p.left !== nothing
        println(io, "left: $(p.left.label)")
    else 
        println(io, "left: nothing")
    end
    if p.right !== nothing
        println(io, "right: $(p.right.label)")
    else
        println(io, "right: nothing")
    end
    println(io, "")
    if p.left !== nothing
        println(io, p.left)
    end
    if p.right !== nothing
        println(io, p.right)
    end
end

function generate_tree(code)
    used_vec = Vector{typeof(code.eins.ixs[1][1])}()
    _generate_tree(code, used_vec)
end

function _generate_tree(code, used_vec)
    # used_vec = Vector{typeof(code.eins.ixs[1])}()
    if !isdefined(code,:eins)
        return nothing
    end
    input_inds = reduce(∪,code.eins.ixs)
    out_inds = code.eins.iy

    inds = setdiff(setdiff(input_inds, used_vec), out_inds)
    if !isempty(inds)
        label = pop!(inds)
        push!(used_vec, label)
        if !isempty(inds)
            left = _generate_tree(code, used_vec)
            if !isempty(inds)
                right = _generate_tree(code, used_vec)
            else
                right = nothing
            end
        else
            used_vec = Vector{typeof(code.eins.ixs[1][1])}()
            left = _generate_tree(code.args[1], used_vec)
            right = _generate_tree(code.args[2], used_vec)
        end
       
        tensor = ones(2,2,2)
        return TensorBinaryTree(tensor, left, right, label)
    end
end
