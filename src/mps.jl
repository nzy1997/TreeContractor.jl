mutable struct LabeledMPS{T<:Number, AT<:AbstractArray{T,3}, LT}
    # rank-3 tensors: left virtual bond, physical bond, right virtual bond
    # left/right virtual bond set to be of dimension 1 at the most left/right
    tensors::Vector{AT}
    labels::Vector{LT}
    label_to_index::Dict{LT, Int}
    center::Int
    function LabeledMPS(tensors::Vector{AT},labels::Vector{LT}) where {T<:Number, AT<:AbstractArray{T,3}, LT}
        @assert length(tensors) == length(labels) "LabeledMPS must have the same number of tensors and labels"
        @assert size(tensors[1], 1) == 1 "Left virtual bond must have dimension 1"
        @assert size(tensors[end], 3) == 1 "Right virtual bond must have dimension 1"
        label_to_index = Dict{LT, Int}(zip(labels, 1:length(labels)))
        new{T, AT, LT}(tensors, labels, label_to_index, -1)
    end
end

nsite(mps::LabeledMPS) = length(mps.labels)
num_of_elements(mps::LabeledMPS) = sum(length, mps.tensors)

function code2mps(code::DynamicNestedEinsum{LT}, size_dict::Dict{LT, Int}) where LT
    # labels = Vector{LT}()
    labels = copy(code.eins.iy)
    apply_vec = Vector{Int}()
    tensor_labels = Vector{Vector{LT}}()
    vanish_labels_vec = Vector{Vector{LT}}()
    _code2mps!(code, labels, apply_vec,tensor_labels, vanish_labels_vec)
    tensors = [ones(Float64,1,size_dict[l],1) for l in labels]
    @assert length(tensors) == length(labels) "tensors and labels must have the same length"
    return LabeledMPS(tensors, labels), apply_vec, tensor_labels, vanish_labels_vec
end

function _code2mps!(code::DynamicNestedEinsum{LT}, labels::Vector{LT}, apply_vec::Vector{Int}, tensor_labels::Vector{Vector{LT}}, vanish_labels_vec::Vector{Vector{LT}}) where LT
    @assert !OMEinsum.isleaf(code) "code is a empty code"
    input_inds = reduce(∪,code.eins.ixs)
    out_inds = code.eins.iy
    vanish_labels = setdiff(input_inds, out_inds)

    if !OMEinsum.isleaf(code.args[1])
        _code2mps!(code.args[1], labels, apply_vec, tensor_labels, vanish_labels_vec)
    else
        push!(apply_vec, code.args[1].tensorindex)
        push!(tensor_labels, code.eins.ixs[1])
        push!(vanish_labels_vec, setdiff(setdiff(code.eins.ixs[1], out_inds),code.eins.ixs[2]))
    end

    for label in vanish_labels
        push!(labels, label)
    end
    @assert OMEinsum.isleaf(code.args[2]) "code.args[2] is not a leaf"
    push!(apply_vec, code.args[2].tensorindex)
    push!(tensor_labels, code.eins.ixs[2])
    push!(vanish_labels_vec, vanish_labels)
    return nothing
end

function apply_tensors!(mps::LabeledMPS, apply_vec::Vector{Int}, tensors::Vector{<:AbstractArray{T}}, tensor_labels::Vector{Vector{LT}}, vanish_labels_vec::Vector{Vector{LT}}; maxdim = Inf) where {T<:Number, LT}
    for (i,label, vanish_labels) in zip(apply_vec, tensor_labels, vanish_labels_vec)
        mps = apply_tensor!(mps, tensors[i], label, vanish_labels)
        # @info  mps.tensors .|> size
        if maximum(size.(mps.tensors,1)) > maxdim
            compress!(FullCompress(), mps; maxdim)
            # @info "compress"
            # @info  mps.tensors .|> size
        end
    end
    return mps
end

function apply_tensor!(mps::LabeledMPS, tensor::AbstractArray{T}, tensor_label::Vector{LT}, vanish_labels::Vector{LT}) where {T<:Number, LT}
    @assert length(tensor_label) == length(size(tensor)) "tensor_label and the dimension of tensor are not compatible"
    for (i,l) in enumerate(tensor_label)
        @assert size(mps.tensors[mps.label_to_index[l]],2) == size(tensor,i) "tensor_label and the dimension of tensor are not compatible at index $i"
    end
    
    sorted_tensor_label = sort(1:length(tensor_label); by= x -> mps.label_to_index[tensor_label[x]])
    mps_vec, bd_vec = tensor2mps(permutedims(tensor,sorted_tensor_label))
   
    nsite_mps = nsite(mps)
    pos = 1
    vanish_pos = Int[]
    for i in mps.label_to_index[tensor_label[sorted_tensor_label[1]]]:mps.label_to_index[tensor_label[sorted_tensor_label[end]]]
        if mps.labels[i] == tensor_label[sorted_tensor_label[pos]]
            merge_tensor = mps_vec[pos]
            pos += 1
            if mps.labels[i] ∉ vanish_labels
                mps.tensors[i] = apply_rank_3_tensor(mps.tensors[i], merge_tensor)
            else
                push!(vanish_pos, i)
                mps.tensors[i] = apply_rank_3_tensor_with_vanish(mps.tensors[i], merge_tensor)
                
            end
        else
            mps.tensors[i] = apply_rank_3_tensor(mps.tensors[i], delta_mps(bd_vec[pos], size(mps.tensors[i],2), T))
        end

        # if i > 1
        #     @assert size(mps.tensors[i],1) == size(mps.tensors[i-1],3)
        # end
    end

    if !isempty(vanish_pos)
        for i in vanish_pos
            if i < nsite_mps
                mps.tensors[i+1] = ein"ab,bcd->acd"(mps.tensors[i][:,:,1], mps.tensors[i+1])
            else
                most_right_tensor = findfirst(x -> mps.labels[x] ∉ vanish_labels, nsite_mps:-1:1)
                if isnothing(most_right_tensor)
                    return LabeledMPS([mps.tensors[end]], [-1])
                else
                    most_right_tensor = nsite_mps - most_right_tensor + 1
                    mps.tensors[most_right_tensor] = ein"abc,cd->abd"(mps.tensors[most_right_tensor], mps.tensors[i][:,:,1])
                end
            end
        end

        deleteat!(mps.tensors, vanish_pos)
        deleteat!(mps.labels, vanish_pos)
        mps.label_to_index = Dict(zip(mps.labels, 1:length(mps.labels)))
        mps.center = -1
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

#          o-g-        
#   d| e/ f|     =>   d|       g/
# -a-o--b--o-c-     -a-o--be--o-c-
function apply_rank_3_tensor_with_vanish(tensor::AbstractArray{T,3}, merge_tensor::AbstractArray{T,3}) where T
    merge_code = ein"bfc,efg->becg"
    m = merge_code(tensor, merge_tensor)
    return reshape(m, size(m,1)*size(m,2), size(m,3)*size(m,4),1)
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


function contract_with_mps(optcode::DynamicNestedEinsum{LT}, tensors::Vector{<:AbstractArray{T}}, size_dict::Dict{LT, Int};maxdim = Inf) where {T<:Number, LT}
    mps, apply_vec, tensor_labels, vanish_labels_vec = code2mps(optcode, size_dict)
    mps = apply_tensors!(mps, apply_vec, tensors, tensor_labels, vanish_labels_vec; maxdim)
    return mps.tensors
end

function random_mps(::Type{T}, N::Int; maxdim::Int, d::Int=2, amplitude::Real=1.0) where T
    @assert N > 0 "Number of sites must be positive, got: $N"
    @assert maxdim > 0 "Maximum bond dimension must be positive, got: $maxdim"
    @assert d > 0 "Physical dimension must be greater than 0, got: $d"
    return LabeledMPS([T(amplitude) .* randn(T, min(d^(i-1), d^(N-i+1), maxdim), d, min(d^i, d^(N-i), maxdim)) for i in 1:N], [i for i in 1:N])
end

function Base.vec(mps::LabeledMPS{T,AT,LT}) where {T,AT,LT<:Integer}
    @assert nsite(mps) <= 30 "MPS too large to be converted to a vector"

    max_label = maximum(mps.labels)
    ixs = Vector{Int}[]
    iy = Int[]
    max_label +=1
    previdx = max_label
    for k in 1:nsite(mps)
        physical = mps.labels[k]
        max_label +=1
        nextidx = max_label
        push!(ixs, [previdx, physical, nextidx])
        k == 1 && push!(iy, previdx)
        push!(iy, physical)
        k == nsite(mps) && push!(iy, nextidx)
        previdx = nextidx
    end
    size_dict = OMEinsum.get_size_dict(ixs, mps.tensors)
    code = optimize_code(
        DynamicEinCode(ixs, iy), size_dict, GreedyMethod()
    )

    return vec(code(mps.tensors...))
end