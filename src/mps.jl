struct LabeledMPS{T<:Number, AT<:AbstractArray{T,3}, LT}
    # rank-3 tensors: left virtual bond, physical bond, right virtual bond
    # left/right virtual bond set to be of dimension 1 at the most left/right
    tensors::Vector{AT}
    labels::Vector{LT}
    label_to_index::Dict{LT, Int}
    function LabeledMPS(tensors::Vector{AT},labels::Vector{LT}) where {T<:Number, AT<:AbstractArray{T,3}, LT}
        @assert length(tensors) == length(labels) "LabeledMPS must have the same number of tensors and labels"
        @assert length(tensors) > 1 "LabeledMPS must have at least 2 tensors"
        @assert size(tensors[1], 1) == 1 "Left virtual bond must have dimension 1"
        @assert size(tensors[end], 3) == 1 "Right virtual bond must have dimension 1"
        label_to_index = Dict{LT, Int}(zip(labels, 1:length(labels)))
        new{T, AT, LT}(tensors, labels, label_to_index)
    end
end

function code2mps(code::DynamicNestedEinsum{LT}, size_dict::Dict{LT, Int}) where LT
    @assert isempty(code.eins.iy) "Output labels must be empty"
    labels = Vector{LT}()
    apply_vec = Vector{Int}()
    tensor_labels = Vector{Vector{LT}}()
    _code2mps!(code, labels, apply_vec,tensor_labels)
    tensors = [ones(Float64,1,size_dict[l],1) for l in labels]
    @assert length(tensors) == length(labels) "tensors and labels must have the same length"
    return LabeledMPS(tensors, labels), apply_vec, tensor_labels
end

function _code2mps!(code::DynamicNestedEinsum{LT}, labels::Vector{LT}, apply_vec::Vector{Int}, tensor_labels::Vector{Vector{LT}}) where LT
    if !isdefined(code,:eins)
        error("code is a empty code")
        push!(apply_vec, code.tensorindex)
        return nothing
    end
    input_inds = reduce(∪,code.eins.ixs)
    out_inds = code.eins.iy

    for label in input_inds
        if label ∉ out_inds
            push!(labels, label)
        end
    end
    push!(apply_vec, code.args[2].tensorindex)
    push!(tensor_labels, code.eins.ixs[2])
    if isdefined(code.args[1],:eins)
        _code2mps!(code.args[1], labels, apply_vec, tensor_labels)
        
    else
        push!(apply_vec, code.args[1].tensorindex)
        push!(tensor_labels, code.eins.ixs[1])
    end

    return nothing
end

function apply_tensors!(mps::LabeledMPS, apply_vec::Vector{Int}, tensors::Vector{<:AbstractArray{T}}, tensor_labels::Vector{Vector{LT}}) where {T<:Number, LT}
    for i in length(apply_vec):-1:1
        apply_tensor!(mps, tensors[i], tensor_labels[i])
        @show  mps.tensors .|> size
        # return
    end
    return mps
end

function apply_tensor!(mps::LabeledMPS, tensor::AbstractArray{T}, tensor_label::Vector{LT}) where {T<:Number, LT}
    @assert length(tensor_label) == length(size(tensor)) "tensor_label and the dimension of tensor are not compatible"
    
    sorted_tensor_label = sort(1:length(tensor_label); by= x -> mps.label_to_index[tensor_label[x]])

    mps_vec, bd_vec = tensor2mps(permutedims(tensor,sorted_tensor_label))
   
    
    pos = 1
    for i in mps.label_to_index[tensor_label[sorted_tensor_label[1]]]:mps.label_to_index[tensor_label[sorted_tensor_label[end]]]
        if mps.labels[i] == tensor_label[sorted_tensor_label[pos]]
            merge_tensor = mps_vec[pos]
            pos += 1
        else
            merge_tensor = delta_mps(bd_vec[pos], size(mps.tensors[i],2), T)
        end
        mps.tensors[i] = apply_rank_3_tensor(mps.tensors[i], merge_tensor)

        # if i > 1
            # @assert size(mps.tensors[i],1) == size(mps.tensors[i-1],3)
        # end
    end
     return mps
end

#          o-g-        
#   d| e/ f|     =>   d|     f| g/
# -a-o--b--o-c-     -a-o--be--o-c-
function apply_rank_3_tensor(tensor::AbstractArray{T,3}, merge_tensor::AbstractArray{T,3}) where T
    merge_code = ein"bfc,efg->befcg"
    m = merge_code(tensor, merge_tensor)
    return reshape(m, size(m,1)*size(m,2), size(m,3), size(m,4)*size(m,5))
end

function tensor2mps(tensor::AbstractArray{T}) where T
    tensor_size = tensor |> size
    mps_vec = Vector{Array{T,3}}()
    bd_vec = Vector{Int}(undef, length(tensor_size)+1)
    bd_vec[1] = 1
    for (i,s) in enumerate(tensor_size[1:end-1])
        mat = reshape(tensor, s*bd_vec[i], :)
        U,S,V = svd(mat)
        tensor = Diagonal(S)*V'
        bd_vec[i+1] = size(U,2)
        push!(mps_vec, reshape(U, bd_vec[i], s, bd_vec[i+1]))
    end
    push!(mps_vec, reshape(tensor, bd_vec[end-1], tensor_size[end], 1))
    bd_vec[end] = 1
    return mps_vec, bd_vec
end


function delta_mps(bd::Int, n::Int,T::Type{<:Number})
    tn = zeros(T,bd,n,bd)
    for i in 1:n
        tn[:,i,:] = I(bd)
    end
    return tn
end

#   d|      =>   
# -a-o-b-        -a-o-b-


#         e|     =>   
# -a-o--b--o-c-      a-o-c-
function contract_mps(tensors::Vector{<:AbstractArray{T,3}}) where T
    tensor = ein"adb->ab"(tensors[1])
    contract_code = ein"ab,bec->ac"
    for t in tensors[2:end]
        tensor = contract_code(tensor, t)
    end
    return tensor
end
contract_mps(mps::LabeledMPS) = contract_mps(mps.tensors)


function contract_with_mps(optcode::DynamicNestedEinsum{LT}, tensors::Vector{<:AbstractArray{T}}, size_dict::Dict{LT, Int}) where {T<:Number, LT}
    mps, apply_vec, tensor_labels = code2mps(optcode, size_dict)
    mps = apply_tensors!(mps, apply_vec, tensors, tensor_labels)
    return contract_mps(mps.tensors)[]
end