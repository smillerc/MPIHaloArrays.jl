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
 - `com_model`: [p2p] Communication model, e.g. :p2p is point-to-point (Isend, Irecv), :rma is onesided (Get,Put), :shared is MPI's shared memory model
"""
function MPIHaloArray(A::AbstractArray{T,N}, topo::CartesianTopology, nhalo::Int; do_corners=true, com_model=:p2p) where {T,N}
    local_di = Vector{DataIndices}(undef, N)
    global_di = Vector{DataIndices}(undef, N)

    A_with_halo = pad_with_halo(A, nhalo)

    if topo.dimension > length(size(A))
        @error "Dimensionality of the AbstractParallelTopology ($(topo.dimension)) > the dimensionality of the array A ($(length(size(A))))"
    end

    for dim in 1:N
        lo_halo_start, lo_halo_end, lo_dom_start, lo_dom_end = lo_indices(A_with_halo, dim, nhalo)
        hi_dom_start, hi_dom_end, hi_halo_start, hi_halo_end = hi_indices(A_with_halo, dim, nhalo)

        local_di[dim] = DataIndices((lo_halo_start, lo_halo_end), 
                                    (lo_dom_start , lo_dom_end),
                                    (lo_dom_start , hi_dom_end),
                                    (hi_dom_start , hi_dom_end),
                                    (hi_halo_start, hi_halo_end))
        
        # TODO: this is wrong -- it assumes a constant size
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

function filldomain!(A::MPIHaloArray, fillval)
    domain = domainview(A)
    fill!(domain, fillval)
end

"""Return a `view` of the domain data within the `A::MPIHaloArray`"""
function domainview(A::MPIHaloArray)
    li = localindices(A)
    lo = li[1:2:end] # low indicies
    hi = li[2:2:end] # high indices
    viewranges = [UnitRange(l,h) for (l,h) in zip(lo, hi)] |> Tuple
    view(A.data, viewranges...)
end

function localindices(A::MPIHaloArray{T,1}) where {T}
    ilo, ihi = A.local_indices[1].domain
    return (ilo, ihi)
end

function localindices(A::MPIHaloArray{T,2}) where {T}
    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain
    return (ilo, ihi, jlo, jhi)
end

function localindices(A::MPIHaloArray{T,3}) where {T}
    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain
    klo, khi = A.local_indices[3].domain
    return (ilo, ihi, jlo, jhi, klo, khi)
end

function globalindices(A::MPIHaloArray{T,1}) where {T}
    ilo, ihi = A.global_indices[1].domain
    return (ilo, ihi)
end

function globalindices(A::MPIHaloArray{T,2}) where {T}
    ilo, ihi = A.global_indices[1].domain
    jlo, jhi = A.global_indices[2].domain
    return (ilo, ihi, jlo, jhi)
end

function globalindices(A::MPIHaloArray{T,3}) where {T}
    ilo, ihi = A.global_indices[1].domain
    jlo, jhi = A.global_indices[2].domain
    klo, khi = A.global_indices[3].domain
    return (ilo, ihi, jlo, jhi, klo, khi)
end

@inline function lo_indices(A::MPIHaloArray{T,N}, dim) where {T,N}
    lo_indices(A.data, dim, A.nhalo)
end

@inline function hi_indices(A::MPIHaloArray{T,N}, dim) where {T,N}
    hi_indices(A.data, dim, A.nhalo)
end

end
