module MPIHaloArrays

using MPI
using OffsetArrays

include("topology.jl")
include("partitioning.jl")

# using .ParallelTopologies

export MPIHaloArray
export ParallelTopology, CartesianTopology
export neighbor, neighbors
export ilo_neighbor, ihi_neighbor, jlo_neighbor, jhi_neighbor, klo_neighbor, khi_neighbor
export lo_indices, hi_indices, fillhalo!, filldomain!
export sync_edges!
export scatterglobal, gatherglobal

struct DataIndices{T <: Integer}
    lo_halo::NTuple{2,T}               # (start,end) indices of the halo region on the low side
    lo_halo_domain_donor::NTuple{2,T}  # (start,end) indices of the domain region on the low side that "donates" to the halo region
    domain::NTuple{2,T}                # (start,end) indices of the real domain
    hi_halo_domain_donor::NTuple{2,T}  # (start,end) indices of the domain region on the high side that "donates" to the halo region
    hi_halo::NTuple{2,T}               # (start,end) indices of the halo region on the high side
end

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
mutable struct MPIHaloArray{T,N} <: AbstractArray{T,N}
    data::AbstractArray{T,N}
    nhalo::Int
    rank::Int
    topology::CartesianTopology
    local_indices::Vector{DataIndices}
    global_indices::Vector{DataIndices}
    do_corners::Bool
    # MPIHaloArray{T}(sizes::Vararg{<:Integer,N}, nhalo) where {T,N} = MPIHaloArray(Array{T, 2}(undef, sizes...), nhalo, 0...)
end

include("utils/indexing.jl")
include("sync_edges.jl")
include("scattergather.jl")

"""MPIHaloArray constructor

# Arguments
 - `A`: AbstractArray{T,N}
 - `topo`: Parallel topology type, e.g. CartesianTopology
 - `nhalo`: Number of halo cells

# Keyword Arguments
 - `do_corners`: [true] Exchange corner halo regions 
 - `com_model`: [p2p] Communication model, e.g. :p2p is point-to-point (Isend, Irecv), :rma is onesided (Get,Put), :shared is MPI's shared memory model
"""
function MPIHaloArray(A::AbstractArray{T,N}, topo::CartesianTopology, nhalo::Int; do_corners=true, com_model=:p2p) where {T,N}
    local_di = Vector{DataIndices}(undef, N)
    global_di = Vector{DataIndices}(undef, N)

    A_with_halo = pad_with_halo(A, nhalo)

    if length(topo.global_dims) > length(size(A))
        @error "Dimensionality of the ParallelTopology ($(length(topo.global_dims))) > the dimensionality of the array A ($(length(size(A))))"
    end

    for dim in 1:N
        lo_halo_start, lo_halo_end, lo_dom_start, lo_dom_end = lo_indices(A_with_halo, dim, nhalo)
        hi_dom_start, hi_dom_end, hi_halo_start, hi_halo_end = hi_indices(A_with_halo, dim, nhalo)

        local_di[dim] = DataIndices((lo_halo_start, lo_halo_end), 
                                    (lo_dom_start , lo_dom_end),
                                    (lo_dom_start , hi_dom_end),
                                    (hi_dom_start , hi_dom_end),
                                    (hi_halo_start, hi_halo_end))
        
        global_offset = (hi_halo_end - lo_halo_start + 1) * topo.coords[dim]

        global_di[dim] = DataIndices((lo_halo_start, lo_halo_end) .+ global_offset, 
                                     (lo_dom_start , lo_dom_end ) .+ global_offset,
                                     (lo_dom_start , hi_dom_end ) .+ global_offset,
                                     (hi_dom_start , hi_dom_end ) .+ global_offset,
                                     (hi_halo_start, hi_halo_end) .+ global_offset)
    end

    update_halo_data!(A, A_with_halo, local_di)
    MPIHaloArray(A_with_halo, nhalo, topo.rank, topo, local_di, global_di, do_corners)
end

function update_halo_data!(A_no_halo::AbstractArray{T,1}, A_with_halo::AbstractArray{T,1}, local_data_indices) where {T}
    ilo_dom, ihi_dom = local_data_indices[1].domain
    ilo, ihi = firstindex(A_no_halo), lastindex(A_no_halo)
    for (i_nh, i_wh) in zip(ilo:ihi, ilo_dom:ihi_dom)
        A_with_halo[i_wh] = A_no_halo[i_nh]
    end
end

function update_halo_data!(A_no_halo::AbstractArray{T,2}, A_with_halo::AbstractArray{T,2}, local_data_indices) where {T}
    ilo_dom, ihi_dom = local_data_indices[1].domain
    jlo_dom, jhi_dom = local_data_indices[2].domain

    CI = CartesianIndices(A_no_halo)
    ilo, jlo = first(CI) |> Tuple
    ihi, jhi = last(CI) |> Tuple

    for (j_nh, j_wh) in zip(jlo:jhi, jlo_dom:jhi_dom)
        for (i_nh, i_wh) in zip(ilo:ihi, ilo_dom:ihi_dom)
            A_with_halo[i_wh, j_wh] = A_no_halo[i_nh, j_nh]
        end
    end
end

function update_halo_data!(A_no_halo::AbstractArray{T,3}, A_with_halo::AbstractArray{T,3}, local_data_indices) where {T}
    ilo_dom, ihi_dom = local_data_indices[1].domain
    jlo_dom, jhi_dom = local_data_indices[2].domain
    klo_dom, khi_dom = local_data_indices[3].domain

    CI = CartesianIndices(A_no_halo)
    ilo, jlo, klo = first(CI) |> Tuple
    ihi, jhi, khi = last(CI) |> Tuple

    for (k_nh, k_wh) in zip(klo:khi, klo_dom:khi_dom)
        for (j_nh, j_wh) in zip(jlo:jhi, jlo_dom:jhi_dom)
            for (i_nh, i_wh) in zip(ilo:ihi, ilo_dom:ihi_dom)
                A_with_halo[i_wh, j_wh, k_wh] = A_no_halo[i_nh, j_nh, k_nh]
            end
        end
    end
end

function pad_with_halo(A,nhalo)
    new_dims = size(A) .+ 2nhalo
    A_new = similar(A, new_dims)
    fill!(A_new, zero(eltype(A)))
    return A_new
end

# Required interface overloads to be an AbstractArray
Base.size(A::MPIHaloArray) = size(A.data)
Base.getindex(A::MPIHaloArray{T,N}, i::Int) where {T,N} = getindex(A.data, i)
Base.getindex(A::MPIHaloArray{T,N}, I::Vararg{Int, N}) where {T,N} = getindex(A.data, I...)
Base.setindex!(A::MPIHaloArray{T,N}, v, i::Int) where {T,N} = setindex!(A.data, v, i)
Base.setindex!(A::MPIHaloArray{T,N}, v, I::Vararg{Int, N}) where {T,N} = setindex!(A.data, v..., I...)

function fillhalo!(A::MPIHaloArray{T,1}, fillval) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo

    iloA = @view A.data[ilo_halo_start:ilo_halo_end, :]
    ihiA = @view A.data[ihi_halo_start:ihi_halo_end, :]

    fill!(iloA, fillval)
    fill!(ihiA, fillval)
end

function fillhalo!(A::MPIHaloArray{T,2}, fillval) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo

    iloA = @view A.data[ilo_halo_start:ilo_halo_end, :]
    ihiA = @view A.data[ihi_halo_start:ihi_halo_end, :]
    jloA = @view A.data[:, jlo_halo_start:jlo_halo_end]
    jhiA = @view A.data[:, jhi_halo_start:jhi_halo_end]

    fill!(iloA, fillval)
    fill!(ihiA, fillval)
    fill!(jloA, fillval)
    fill!(jhiA, fillval)
end

function fillhalo!(A::MPIHaloArray{T,3}, fillval) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    klo_halo_start, klo_halo_end = A.local_indices[3].lo_halo
    khi_halo_start, khi_halo_end = A.local_indices[3].hi_halo

    iloA = @view A.data[ilo_halo_start:ilo_halo_end, :, :]
    ihiA = @view A.data[ihi_halo_start:ihi_halo_end, :, :]
    jloA = @view A.data[:, jlo_halo_start:jlo_halo_end, :]
    jhiA = @view A.data[:, jhi_halo_start:jhi_halo_end, :]
    kloA = @view A.data[:, :, klo_halo_start:klo_halo_end]
    khiA = @view A.data[:, :, khi_halo_start:khi_halo_end]

    fill!(iloA, fillval)
    fill!(ihiA, fillval)
    fill!(jloA, fillval)
    fill!(jhiA, fillval)
    fill!(kloA, fillval)
    fill!(khiA, fillval)
end

function filldomain!(A::MPIHaloArray{T,1}, fillval) where {T}
    ilo, ihi  = A.local_indices[1].domain
    domA = @view A.data[ilo:ihi]
    fill!(domA, fillval)
end

function filldomain!(A::MPIHaloArray{T,2}, fillval) where {T}
    ilo, ihi  = A.local_indices[1].domain
    jlo, jhi  = A.local_indices[2].domain
    domA = @view A.data[ilo:ihi, jlo:jhi]
    fill!(domA, fillval)
end

function filldomain!(A::MPIHaloArray{T,3}, fillval) where {T}
    ilo, ihi  = A.local_indices[1].domain
    jlo, jhi  = A.local_indices[2].domain
    klo, khi  = A.local_indices[3].domain
    domA = @view A.data[ilo:ihi, jlo:jhi, klo:khi]
    fill!(domA, fillval)
end


# include("ops.jl")

# """
#     applylocalfunc(f, A::MPIHaloArray)
# Execute the function f on the part of a owned by the current rank. It is assumed f does not modify the local part.
# """
# function applylocalfunc(f, A::MPIHaloArray)
#     MPI.Win_lock(MPI.LOCK_SHARED, A.rank, 0, A.window)
#     result = f(A.data)
#     MPI.Win_unlock(A.rank, A.window)
#     return result
# end

end
