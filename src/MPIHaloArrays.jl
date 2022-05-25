module MPIHaloArrays

using MPI
using OffsetArrays

include("topology.jl")
include("partitioning.jl")
include("utils/dataindices.jl")

export MPIHaloArray
export AbstractParallelTopology, CartesianTopology
export neighbor, neighbors
export ilo_neighbor, ihi_neighbor, jlo_neighbor, jhi_neighbor, klo_neighbor, khi_neighbor
export lo_indices, hi_indices, fillhalo!, filldomain!
export updatehalo!, domainview
export scatterglobal, gatherglobal
export localindices, globalindices
export globalmin, globalmax, globalsum

"""
MPIHaloArray

# Fields
 - `data`: AbstractArray{T,N} - contains the local data on the current rank
 - `partitioning`: partitioning datatype
 - `comm`: MPI communicator
 - `window`: MPI window
 - `neighbor_ranks` : Vector{Int} - IDs of the neighboring arrays/MPI procs
 - `coords` : Vector{Int} - Coordinates in the global MPI space
 - `rank`: Current MPI rank
"""
mutable struct MPIHaloArray{T,N,NN} <: AbstractArray{T,N}
    data::AbstractArray{T,N}
    nhalo::Int
    rank::Int
    topology::CartesianTopology
    local_indices::Vector{DataIndices{Int}}
    global_indices::Vector{DataIndices{Int}}
    do_corners::Bool
    halo_dims::NTuple{NN,Int64}
end


include("utils/indexing.jl")
include("updatehalo.jl")
include("scattergather.jl")
include("ops.jl")

"""MPIHaloArray constructor

# Arguments
 - `A`: AbstractArray{T,N}
 - `topo`: Parallel topology type, e.g. CartesianTopology
 - `nhalo`: Number of halo cells

# Keyword Arguments
 - `do_corners`: [true] Exchange corner halo regions
 - `com_model`: [:p2p] Communication model, e.g. :p2p is point-to-point (Isend, Irecv), :rma is onesided (Get,Put), :shared is MPI's shared memory model
"""
function MPIHaloArray(A::AbstractArray{T,NN}, topo::CartesianTopology, nhalo::Int, halo_dims::NTuple{N,Int}; do_corners=true, com_model=:p2p) where {T,N,NN}

    if any(halo_dims .> ndims(A))
        error("Some (or all) of the given halo_dims $(halo_dims) are incompatible with the dimensionality of the given array A")
    end

    if topo.dimension != length(halo_dims)
        error("Mismatched topology/halo exchange dimensionality (topology is $(topo.dimension)D, halo exchange dimensions are $(length(halo_dims))D)")
    end

    local_di = Vector{DataIndices{Int64}}(undef, NN)
    global_di = Vector{DataIndices{Int64}}(undef, NN)

    A_with_halo = pad_with_halo(A, nhalo, halo_dims)

    for dim in 1:ndims(A)

        if dim in halo_dims
            lo_halo_start, lo_halo_end, lo_dom_start, lo_dom_end = lo_indices(A_with_halo, dim, nhalo)
            hi_dom_start, hi_dom_end, hi_halo_start, hi_halo_end = hi_indices(A_with_halo, dim, nhalo)

            local_di[dim] = DataIndices((lo_halo_start, lo_halo_end),
                (lo_dom_start, lo_dom_end),
                (lo_dom_start, hi_dom_end),
                (hi_dom_start, hi_dom_end),
                (hi_halo_start, hi_halo_end))

            # TODO: this is wrong -- it assumes a constant size
            global_offset = (hi_halo_end - lo_halo_start + 1) * topo.coords[dim]
            global_di[dim] = DataIndices((lo_halo_start, lo_halo_end) .+ global_offset,
                (lo_dom_start, lo_dom_end) .+ global_offset,
                (lo_dom_start, hi_dom_end) .+ global_offset,
                (hi_dom_start, hi_dom_end) .+ global_offset,
                (hi_halo_start, hi_halo_end) .+ global_offset)
        else
            lo = firstindex(A, dim)
            hi = lastindex(A, dim)
            local_di[dim] = DataIndices((lo, lo),
                (lo, lo),
                (lo, hi),
                (hi, hi),
                (hi, hi))

            global_offset = (hi - lo + 1) * topo.coords[dim]
            global_di[dim] = DataIndices((lo, lo) .+ global_offset,
                (lo, lo) .+ global_offset,
                (lo, hi) .+ global_offset,
                (hi, hi) .+ global_offset,
                (hi, hi) .+ global_offset)
        end
    end

    update_halo_data!(A, A_with_halo, halo_dims, nhalo)
    MPIHaloArray(A_with_halo, nhalo, topo.rank, topo, local_di, global_di, do_corners, halo_dims)
end


function MPIHaloArray(A::AbstractArray{T,N}, topo::CartesianTopology, nhalo::Int, halo_dims::Int; do_corners=true, com_model=:p2p) where {T,N}
    MPIHaloArray(A, topo, nhalo, tuple(halo_dims); do_corners=do_corners, com_model=com_model)
end

function MPIHaloArray(A::AbstractArray{T,N}, topo::CartesianTopology, nhalo::Int; do_corners=true, com_model=:p2p) where {T,N}
    halo_dims = Tuple(1:ndims(A)) # unless specified, assume that halo exchange is done on every dimension
    MPIHaloArray(A, topo, nhalo, halo_dims; do_corners=do_corners, com_model=com_model)
end


"""
    update_halo_data!(A_no_halo, A_with_halo, halo_dims, nhalo)

During the construction of an `MPIHaloArray`, the data must be padded by the number of halo regions in
each respective dimension. This function copies the data from `A_no_halo` (which is the original array)
to `A_with_halo` (which is the underlying array within the `MPIHaloArray`)

# Arguments
 - `A_no_halo::Array`: Array without the halo regions
 - `A_with_halo::Array`: Array padded with the halo regions
 - `halo_dims::Tuple`: Dimensions that halo exchanges take place on
 - `nhalo::Int`: Number of halo cells in each respective dimension
"""
function update_halo_data!(A_no_halo, A_with_halo, halo_dims, nhalo)

    view_dims = Vector{UnitRange{Int64}}(undef, ndims(A_no_halo))

    # Construct the ranges in each dimension that the real data lives within
    for dim in 1:ndims(A_no_halo)
        if dim in halo_dims
            _, _, lo_dom_start, _ = lo_indices(A_with_halo, dim, nhalo)
            _, hi_dom_end, _, _ = hi_indices(A_with_halo, dim, nhalo)

            view_dims[dim] = lo_dom_start:hi_dom_end
        else
            view_dims[dim] = axes(A_no_halo, dim)
        end
    end

    A_with_halo_data = @view A_with_halo[view_dims...]
    copy!(A_with_halo_data, A_no_halo)

    return nothing
end

"""
	pad_with_halo(A, nhalo, halo_dims)

Increase the size of the array `A` along the halo exchange dimensions to make room for the new halo regions

# Arguments
 - `A::AbstractArray`: Array to increase in size
 - `nhalo::Int`, number of halo cells along each dimension, e.g., 2
 - `halo_dims::Tuple`: Set of dimensions to do halo exchange along
"""
function pad_with_halo(A, nhalo, halo_dims)
    # use a mutable vector since some dims may change size
    new_dims = size(A) |> collect

    for d in halo_dims
        new_dims[d] += 2nhalo
    end

    A_new = similar(A, new_dims |> Tuple)
    fill!(A_new, zero(eltype(A)))

    return A_new
end

# Required interface overloads to be an AbstractArray
Base.size(A::MPIHaloArray) = size(A.data)
Base.getindex(A::MPIHaloArray{T,N}, i::Int) where {T,N} = getindex(A.data, i)
Base.getindex(A::MPIHaloArray{T,N}, I::Vararg{Int,N}) where {T,N} = getindex(A.data, I...)
Base.setindex!(A::MPIHaloArray{T,N}, v, i::Int) where {T,N} = setindex!(A.data, v, i)
Base.setindex!(A::MPIHaloArray{T,N}, v, I::Vararg{Int,N}) where {T,N} = setindex!(A.data, v..., I...)

"""
    fillhalo!(A::MPIHaloArray, fillvalue)

Fill the halo regions with a particular `fillvalue`

# Arguments
 - `A::MPIHaloArray`
 - `fillvalue`: value to fill the halo regions of A with
"""
function fillhalo!(A::MPIHaloArray, fillvalue)

    for dim in A.halo_dims
        # create a vector of index ranges and selectively change
        # the halo dimension to be the hi/lo subarray, so that we can
        # fill them with the fillvalue
        lo_view_dims = axes(A) .|> UnitRange |> collect
        hi_view_dims = axes(A) .|> UnitRange |> collect

        lo_halo_start, lo_halo_end = A.local_indices[dim].lo_halo
        hi_halo_start, hi_halo_end = A.local_indices[dim].hi_halo

        lo_view_dims[dim] = lo_halo_start:lo_halo_end
        hi_view_dims[dim] = hi_halo_start:hi_halo_end

        loA = @view A[lo_view_dims...]
        hiA = @view A[hi_view_dims...]

        fill!(loA, fillvalue)
        fill!(hiA, fillvalue)
    end

    return nothing
end

"""Fill the domain data with a single `filval`"""
function filldomain!(A::MPIHaloArray, fillval)
    domain = domainview(A)
    fill!(domain, fillval)
end

"""Return a view or `SubArray` of the domain data within the `MPIHaloArray`"""
function domainview(A::MPIHaloArray)
    li = localindices(A)
    lo = li[1:2:end] # low indicies
    hi = li[2:2:end] # high indices
    viewranges = [UnitRange(l, h) for (l, h) in zip(lo, hi)] |> Tuple
    view(A.data, viewranges...)
end

function localindices(A::MPIHaloArray)
     getindices(A.local_indices)
end

function globalindices(A::MPIHaloArray)
     getindices(A.global_indices)
end

function getindices(dataindices)
    indices = zeros(Int, 2length(dataindices))

    for dim in 1:length(dataindices)
        lo, hi = dataindices[dim].domain
        indices[2(dim - 1) + 1] = lo
        indices[2(dim - 1) + 2] = hi
    end

    return tuple(indices...)
end

@inline function lo_indices(A::MPIHaloArray, dim)
    lo_indices(A.data, dim, A.nhalo)
end

@inline function hi_indices(A::MPIHaloArray, dim)
    hi_indices(A.data, dim, A.nhalo)
end

end
