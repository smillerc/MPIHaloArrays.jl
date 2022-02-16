module MPIHaloArrays

using MPI
using OffsetArrays

include("./topology.jl")

# using .ParallelTopologies

export MPIHaloArray
export ParallelTopology, CartesianTopology
export neighbor, neighbors
export ilo_neighbor, ihi_neighbor, jlo_neighbor, jhi_neighbor, klo_neighbor, khi_neighbor
export lo_indices, hi_indices, fillhalo!, filldomain!
export sync_edges!

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
 - `data`: Array{T,N} - contains the local data on the current rank
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
    # MPIHaloArray{T}(sizes::Vararg{<:Integer,N}, nhalo) where {T,N} = MPIHaloArray(Array{T, 2}(undef, sizes...), nhalo, 0...)
end

# include("sync_edges.jl")
include("utils/indexing.jl")

function MPIHaloArray(A::AbstractArray{T,N}, topo::CartesianTopology, nhalo::Int) where {T,N}
    local_di = Vector{DataIndices}(undef, N)
    global_di = Vector{DataIndices}(undef, N)

    A_with_halo = pad_with_halo(A, nhalo)

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
    MPIHaloArray(A_with_halo, nhalo, topo.rank, topo, local_di, global_di)
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

function filldomain!(A::MPIHaloArray{T,2}, fillval) where {T}

    ilo, ihi  = A.local_indices[1].domain
    jlo, jhi  = A.local_indices[2].domain
    domA = @view A.data[ilo:ihi, jlo:jhi]
    fill!(domA, fillval)
end

"""Sync the edges of the array `A` with it's neighbors"""
function sync_edges!(A::MPIHaloArray{T,2}) where {T}
    # Get the start/end indices for subarray/view extraction of the halo regions
    # ndims = length(size(A))
    # ilo_halo_start, ilo_halo_end, ilo_dom_start, _ = lo_indices(nhalo)
    # jlo_halo_start, jlo_halo_end, jlo_dom_start, _ = lo_indices(nhalo)
    # _, ihi_dom_end, ihi_halo_start, _ = hi_indices(A, ndims - 1, nhalo)
    # _, jhi_dom_end, jhi_halo_start, _ = hi_indices(A, ndims, nhalo)

    nhalo = A.nhalo
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start , ilo_dom_end  = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start , ihi_dom_end  = A.local_indices[1].hi_halo_domain_donor

    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    jlo_dom_start , jlo_dom_end  = A.local_indices[2].lo_halo_domain_donor
    jhi_dom_start , jhi_dom_end  = A.local_indices[2].hi_halo_domain_donor
    

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_start + nhalo - 1, jlo_dom_start:jhi_dom_end]
    jlo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_start + nhalo - 1]
    ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
    jlo_halo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]

    # Define the positions w/in the window for Get/Put operations
    ilo_edge_pos = LinearIndices(A.data)[ilo_dom_start, jlo_dom_start]
    jlo_edge_pos = LinearIndices(A.data)[ilo_dom_start, jlo_dom_start]
    ihi_halo_pos = LinearIndices(A.data)[ihi_halo_start, jlo_dom_start]
    jhi_halo_pos = LinearIndices(A.data)[ilo_dom_start, jhi_halo_start]

    # Create the MPI subarray buffers for transfer
    ilo_buf = MPI.Buffer(ilo_edge)
    jlo_buf = MPI.Buffer(jlo_edge)
    ilo_halo_buf = MPI.Buffer(ilo_halo_edge)
    jlo_halo_buf = MPI.Buffer(jlo_halo_edge)

    # Calculate the offset to move the buffers
    ilo_to_ihi_halo = ihi_halo_pos - ilo_edge_pos
    jlo_to_jhi_halo = jhi_halo_pos - jlo_edge_pos
    ihi_halo_to_ilo = ihi_halo_pos - ilo_edge_pos
    jhi_halo_to_jlo = jhi_halo_pos - jlo_edge_pos

    # Halo exchange
    ilo_neighbor_proc = ilo_neighbor(A.topology)
    ihi_neighbor_proc = ihi_neighbor(A.topology)
    jlo_neighbor_proc = jlo_neighbor(A.topology)
    jhi_neighbor_proc = jhi_neighbor(A.topology)

    # ilo_neighbor, jlo_neighbor, _ = global_domain().neighbors[1,:]
    # ihi_neighbor, jhi_neighbor, _ = global_domain().neighbors[2,:]

    win = MPI.Win_create(A.data, A.topology.comm)
    MPI.Win_fence(0, win)

    MPI.Get(ilo_halo_buf, ilo_neighbor_proc, ihi_halo_to_ilo, win) # ilo halo update
    MPI.Get(jlo_halo_buf, jlo_neighbor_proc, jhi_halo_to_jlo, win) # jlo halo update
    MPI.Put(jlo_buf,      jlo_neighbor_proc, jlo_to_jhi_halo, win) # jhi halo update
    MPI.Put(ilo_buf,      ilo_neighbor_proc, ilo_to_ihi_halo, win) # ihi halo update
    
    MPI.Win_fence(0, win)

    return nothing
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
