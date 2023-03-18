# module ParallelTopologies

using MPI
using OffsetArrays
using LinearAlgebra: norm

# export AbstractParallelTopology, CartesianTopology
# export ilo_neighbor, ihi_neighbor, jlo_neighbor, jhi_neighbor, klo_neighbor, khi_neighbor, neighbor, neighbors

"""An abstract AbstractParallelTopology type that is extended by either a CartesianTopology or GraphTopology (future)"""
abstract type AbstractParallelTopology end

global const I = 0
global const J = 0
global const K = 0
global const ILO = -1
global const IHI = +1
global const JLO = -1
global const JHI = +1
global const KLO = -1
global const KHI = +1


"""
CartesianTopology

The CartesianTopology type holds neighbor information, current rank, etc.

# Fields
 - `comm`: MPI commicator object
 - `nprocs`: Number of total processors (global)
 - `rank`: Current rank
 - `coords`: Coordinates in the global space, i.e. `(0,1,1)`
 - `global_dims`: Dimensions of the global domain, i.e. `(4,4)` is a 4x4 global domain
 - `isperiodic`: Vector{Bool}; Perodicity of each dimension, i.e. `(false, true, true)` means y and z are periodic
 - `neighbors`: OffsetArray{Int}; Neighbor ranks (including corners), indexed as `[[ilo, center, ihi], i, j, k]`
"""
struct CartesianTopology <: AbstractParallelTopology
    comm::MPI.Comm
    nprocs::Int
    rank::Int
    dimension::Int
    coords::NTuple{3,Int}    # [i, j, k]; coordinates in the toplogy
    global_dims::NTuple{3,Int} # [i, j, k]; number of domains in each direction
    isperiodic::NTuple{3,Bool} # [i, j, k]; is this dimension periodic?
    neighbors::OffsetArray{Int, 3}   # [[ilo, center, ihi], i, j, k]; defaults to -1 if no neighbor
end

Base.size(C::CartesianTopology) = C.global_dims[1:C.dimension]

#    Neighbor index convention (using offset arrays for index simplicity)
#    An index of (0,0,1) will give the neighbor in the k+1 direction, (-1,0,-1) is in the i-1, j, k-1 direction
#
#      /    (k)
#    /
#   +----- (i)
#   |
#   | (j)
#
#
#              +-------+-------+-------+
#            /       /       /       /  |
#          /       /       /       /    |
#        +-------+-------+-------+      |
#      /       /       /       /  |     +
#    /       /       /       /    |   / |
#   +-------+-------+-------+     | /   |
#   |       |       |       |     +     |
#   |       |       |       |   / |     +
#   |       |       |       | /   |   / |
#   +-------+-------+-------+     | /   |
#   |       |       |       |     +     |
#   |       |       |       |   / |     +
#   |       |       |       | /   |   /
#   +-------+-------+-------+     | /
#   |       |       |       |     +
#   |       |       |       |   /
#   |       |       |       | /
#   +-------+-------+-------+



"""
    CartesianTopology(comm::MPI.Comm, dims, periodicity; canreorder = false)

Create a CartesianTopology type that holds neighbor information, current rank, etc.

# Arguments
 - `dims`: Vector or Tuple setting the dimensions of the domain in each direction, e.g. (4,3) means a total of 12 procs, with 4 in x and 3 in y
 - `periodicity`: Vector or Tuple of bools to set if the domain is periodic along a specific dimension

# Example
```julia

# Create a topology of 4x4 with periodic boundaries in both directions
P = CartesianTopology((4,4), (true, true))
```
"""
function CartesianTopology(comm::MPI.Comm, dims::NTuple{N, Int}, periodicity::NTuple{N, Bool}; canreorder = false) where {N}
    dims = collect(dims)
    periodicity = collect(periodicity)

    # comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    # use reverse(dims) b/c MPI convention in (k, j, i), but the expected order is (i, j, k)
    mpi_dims = reverse(dims)
    mpi_periodicity = reverse(periodicity)

    @assert length(dims) == length(periodicity) "You must specify periodicity (true/false) for each dimension"
    @assert prod(dims) == nprocs "The number of subdomains ($(prod(dims))) doesn't match the number of procs ($(nprocs)) provided"

    comm_cart = MPI.Cart_create(comm, mpi_dims, mpi_periodicity .|> Int, canreorder)
    coords = MPI.Cart_coords(comm_cart) |> reverse

    coords_tuple = vec_to_ntuple(coords)
    dims_tuple = vec_to_ntuple(dims)
    periodicity_tuple = vec_to_ntuple(periodicity)
    topo_dim = length(dims)

    neighbors = OffsetArray(Int(MPI.API.MPI_PROC_NULL[]) * ones(Int8,3,3,3), -1:1, -1:1, -1:1)

    # MPI convention is (k, j, i), or (z, y, x) which is annoying
    if topo_dim == 1
        ilo, ihi = MPI.Cart_shift(comm_cart, 0, 1)

        neighbors[:,  0, 0] = [ilo, rank, ihi]

    elseif topo_dim == 2
        jlo, jhi = MPI.Cart_shift(comm_cart, 0, 1)
        ilo, ihi = MPI.Cart_shift(comm_cart, 1, 1)

        jhi_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JHI)
        jlo_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JLO)
        jhi_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JHI)
        jlo_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JLO)

        neighbors[:,  1, 0] = [jhi_ilo, jhi , jhi_ihi]
        neighbors[:,  0, 0] = [ilo    , rank, ihi]
        neighbors[:, -1, 0] = [jlo_ilo, jlo , jlo_ihi]

    elseif topo_dim == 3
        ilo, ihi = MPI.Cart_shift(comm_cart, 2, 1)
        jlo, jhi = MPI.Cart_shift(comm_cart, 1, 1)
        klo, khi = MPI.Cart_shift(comm_cart, 0, 1)

        jhi_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JHI, K)
        jhi_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JHI, K)
        jlo_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JLO, K)
        jlo_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JLO, K)

        klo_jhi     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, I  , JHI, KLO)
        klo_jlo     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, I  , JLO, KLO)
        klo_ilo     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, J  , KLO)
        klo_ihi     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, J  , KLO)
        klo_jhi_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JHI, KLO)
        klo_jlo_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JLO, KLO)
        klo_jhi_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JHI, KLO)
        klo_jlo_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JLO, KLO)

        khi_jhi     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, I  , JHI, KHI)
        khi_jlo     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, I  , JLO, KHI)
        khi_ilo     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, J  , KHI)
        khi_ihi     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, J  , KHI)
        khi_jhi_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JHI, KHI)
        khi_jlo_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JLO, KHI)
        khi_jhi_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JHI, KHI)
        khi_jlo_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JLO, KHI)

        # Center
        neighbors[:,  1, K] = [jhi_ilo, jhi , jhi_ihi]
        neighbors[:,  0, K] = [ilo    , rank, ihi]
        neighbors[:, -1, K] = [jlo_ilo, jlo , jlo_ihi]

        neighbors[:,  1, KHI] = [khi_jhi_ilo, khi_jhi, khi_jhi_ihi]
        neighbors[:,  0, KHI] = [khi_ilo    , khi    , khi_ihi]
        neighbors[:, -1, KHI] = [khi_jlo_ilo, khi_jlo, khi_jlo_ihi]

        neighbors[:,  1, KLO] = [klo_jhi_ilo, klo_jhi, klo_jhi_ihi]
        neighbors[:,  0, KLO] = [klo_ilo    , klo    , klo_ihi]
        neighbors[:, -1, KLO] = [klo_jlo_ilo, klo_jlo, klo_jlo_ihi]

    end

    CartesianTopology(comm_cart, nprocs, rank, topo_dim, coords_tuple, dims_tuple, periodicity_tuple, neighbors)
end

function CartesianTopology(comm::MPI.Comm, dims::Vector{Int}, periodicity::Vector{Bool}; canreorder = false) where {N}
    CartesianTopology(comm, tuple(dims...), tuple(periodicity...); canreorder = canreorder)
end

"""
    CartesianTopology(comm::MPI.Comm, periodicity::Bool; canreorder = false)

Create CartesianTopology only with the vector of boundary periodicity given. This finds the
optimal sub-domain ordering for the user.
"""
function CartesianTopology(comm::MPI.Comm, nprocs::Int, periodicity::Bool; canreorder = false)
    CartesianTopology(comm, nprocs, tuple(periodicity); canreorder = canreorder)
end

"""
    CartesianTopology(comm::MPI.Comm, ::Tuple{Bool}; canreorder = false)

Create CartesianTopology only with the vector of boundary periodicity given. This finds the
optimal sub-domain ordering for the user.
"""
function CartesianTopology(comm::MPI.Comm, nprocs::Int, periodicity::NTuple{N,Bool}; canreorder = false) where {N}
    if length(periodicity) == 3
        dims = num_3d_tiles(nprocs)
    elseif length(periodicity) == 2
        dims = num_2d_tiles(nprocs)
    elseif length(periodicity) == 1
        dims = tuple(nprocs)
    end

    CartesianTopology(comm, dims, periodicity; canreorder = canreorder)
end

function vec_to_ntuple(v::Vector{T}) where {T <: Number}
    if length(v) < 3
        return vcat(v, zeros(T, 3 - length(v))) |> Tuple
    else
        return v |> Tuple
    end
end

function vec_to_ntuple(v::Vector{Bool})
    if length(v) < 3
        return vcat(v, falses(3 - length(v))) |> Tuple
    else
        return v |> Tuple
    end
end

"""Helper function to find rank based on 3D offsets"""
function offset_coord_to_rank(comm, dims, periods, i_offset::Int, j_offset::Int, k_offset::Int)
    # coords = MPI.Cart_coords(comm) .+ (i_offset, j_offset, k_offset)
    coords = MPI.Cart_coords(comm) .+ (k_offset, j_offset, i_offset)
    coord_to_rank(comm, dims, periods, coords)
end

"""Helper function to find rank based on 2D offsets"""
function offset_coord_to_rank(comm, dims, periods, i_offset::Int, j_offset::Int)
    coords = MPI.Cart_coords(comm) .+ (j_offset, i_offset)
    # coords = MPI.Cart_coords(comm) .+ (i_offset, j_offset)
    coord_to_rank(comm, dims, periods, coords)
end

"""Helper function to find rank based on coordinates"""
function coord_to_rank(comm, dims, periods, coords)

    mpi_coords = coords #|> reverse # mpi uses reverse coordinate convention
    isvalid = true

    for i in 1:length(dims)
        if (!periods[i] && (mpi_coords[i] >= dims[i] || mpi_coords[i] < 0))
            isvalid = false
            break
        end
    end

    if isvalid
        rank = MPI.Cart_rank(comm, mpi_coords)
    else
        rank = Int(MPI.API.MPI_PROC_NULL[])
    end
    rank
end

function Base.show(io::IO, p::AbstractParallelTopology)
    println(io, "AbstractParallelTopology")
    println(io, "rank $(p.rank) of $(p.nprocs)")
    println(io, "coords:", p.coords)
    println(io, "global_dims:", p.global_dims)
    println(io, "isperiodic:", p.isperiodic)
    println(io, "ilo neighbor: $(ilo_neighbor(p))")
    println(io, "ihi neighbor: $(ihi_neighbor(p))")
    println(io, "jlo neighbor: $(jlo_neighbor(p))")
    println(io, "jhi neighbor: $(jhi_neighbor(p))")
    println(io, "klo neighbor: $(klo_neighbor(p))")
    println(io, "khi neighbor: $(khi_neighbor(p))")
end

"""Neighbor rank in the i-1 direction"""
ilo_neighbor(p::CartesianTopology) = p.neighbors[-1, 0, 0]

"""Neighbor rank in the i+1 direction"""
ihi_neighbor(p::CartesianTopology) = p.neighbors[ 1, 0, 0]

"""Neighbor rank in the j-1 direction"""
jlo_neighbor(p::CartesianTopology) = p.neighbors[ 0,-1, 0]

"""Neighbor rank in the j+1 direction"""
jhi_neighbor(p::CartesianTopology) = p.neighbors[ 0, 1, 0]

"""Neighbor rank in the k-1 direction"""
klo_neighbor(p::CartesianTopology) = p.neighbors[ 0, 0,-1]

"""Neighbor rank in the k+1 direction"""
khi_neighbor(p::CartesianTopology) = p.neighbors[ 0, 0, 1]

"""
    neighbor(p::CartesianTopology, i_offset::Int, j_offset::Int, k_offset::Int)

Find the neighbor rank based on the offesets in `(i,j,k)`. This follows the traditional
array index convention rather than MPI's version, so an `i_offset=1` will shift up in
the array indexing.

# Arguments
 - `p`: CartesianTopology type
 - `i_offset`: Offset in the `i` direction
 - `j_offset`: Offset in the `j` direction
 - `k_offset`: Offset in the `k` direction

# Example:
```julia
# Makes a 4x4 domain with periodic boundaries in both dimensions
P = CartesianTopology((4,4), (true, true))

# Find the ihi neighbor
ihi = neighbor(P,+1,0,0)

# Find the upper ihi corner neighbor (ihi and jhi side)
ihijhi_corner = neighbor(P,+1,+1,0)
```

"""
function neighbor(p::CartesianTopology, i_offset::Int, j_offset::Int, k_offset::Int)
    p.neighbors[i_offset, j_offset, k_offset]
end

function neighbor(p::CartesianTopology, i_offset::Int, j_offset::Int)
    p.neighbors[i_offset, j_offset, 0]
end

function neighbor(p::CartesianTopology, i_offset::Int)
    p.neighbors[i_offset, 0, 0]
end

function neighbors(p::CartesianTopology)
    p.neighbors
end


"""Return all common denominators of n"""
function denominators(n::Integer)
    denominators = Vector{Int}(undef, 0)
    for i in 1:n
        if mod(n, i) == 0
            push!(denominators, i)
        end
    end
    return denominators
end

"""Returns the optimal number of tiles in (i,j) given total number of tiles n"""
function num_2d_tiles(n)
    # find all common denominators of the total number of images
    denoms = denominators(n)

    # find all combinations of common denominators
    # whose product equals the total number of images
    dim1 = Vector{Int}(undef, 0)
    dim2 = Vector{Int}(undef, 0)
    for j in 1:length(denoms)
        for i in 1:length(denoms)
            if denoms[i] * denoms[j] == n
                push!(dim1, denoms[i])
                push!(dim2, denoms[j])
            end
        end
    end
    # pick the set of common denominators with the minimal norm
    # between two elements -- rectangle closest to a square
    num_2d_tiles = [dim1[1], dim2[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i]] .- sqrt(n))
        n2 = norm(num_2d_tiles .- sqrt(n))
        if n1 < n2
            num_2d_tiles = [dim1[i], dim2[i]]
        end
    end

    return tuple(num_2d_tiles...)
end

"""Returns the optimal number of tiles in (i,j,k) given total number of tiles n"""
function num_3d_tiles(n)
    # find all common denominators of the total number of images
    denoms = denominators(n)

    # find all combinations of common denominators
    # whose product equals the total number of images
    dim1 = Vector{Int}(undef, 0)
    dim2 = Vector{Int}(undef, 0)
    dim3 = Vector{Int}(undef, 0)
    for k in 1:length(denoms)
        for j in 1:length(denoms)
            for i in 1:length(denoms)
                if denoms[i] * denoms[j] * denoms[k] == n
                    push!(dim1, denoms[i])
                    push!(dim2, denoms[j])
                    push!(dim3, denoms[k])
                end
            end
        end
    end
    # pick the set of common denominators with the minimal norm
    # between two elements -- rectangle closest to a square
    num_3d_tiles = [dim1[1], dim2[1], dim3[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i], dim3[i]] .- sqrt(n))
        n2 = norm(num_3d_tiles .- sqrt(n))
        if n1 < n2
            num_3d_tiles = [dim1[i], dim2[i], dim3[i]]
        end
    end

    return tuple(num_3d_tiles...)
end

# """
# Given an input I dimensions of the total computational domain,
# returns an array of start and end indices [ilo,ihi]
# """
# function tile_indices_1d1(dims::Integer, ntiles::Integer, id::Integer)
#     indices = zeros(Int, 2)
#     tile_size = dims ÷ ntiles

#     # start and end indices assuming equal tile sizes
#     indices[1] = (id - 1) * tile_size + 1
#     indices[2] = indices[1] + tile_size - 1

#     # if we have any remainder, distribute it to the tiles at the end
#     offset = ntiles - mod(dims, ntiles)
#     if id > offset
#         indices[1] = indices[1] + id - offset - 1
#         indices[2] = indices[2] + id - offset
#     end
#     return indices
# end

# """
# Given an input (I,J) dimensions of the total computational domain,
# returns an array of start and end indices [ilo,ihi,jlo,jhi]
# """
# function tile_indices_2d(dims, ntiles::Integer, id::Integer)
#     indices = zeros(Int, 4)
#     tiles = num_2d_tiles(ntiles)
#     tiles_ij = tile_id_to_ij(id, ntiles)
#     indices[1:2] = tile_indices_1d(dims[1], tiles[1], tiles_ij[1])
#     indices[3:4] = tile_indices_1d(dims[2], tiles[2], tiles_ij[2])
#     return indices
# end

# """
# Given an input (I,J,K) dimensions of the total computational domain,
# returns an array of start and end indices [ilo,ihi,jlo,jhi,klo,khi]
# """
# function tile_indices_3d(dims, ntiles::Integer, id::Integer)
#     indices = zeros(Int, 6)
#     tiles = num_3d_tiles(ntiles)
#     tiles_ij = tile_id_to_ijk(id, ntiles)
#     indices[1:2] = tile_indices_1d(dims[1], tiles[1], tiles_ij[1])
#     indices[3:4] = tile_indices_1d(dims[2], tiles[2], tiles_ij[2])
#     indices[5:6] = tile_indices_1d(dims[3], tiles[3], tiles_ij[3])
#     return indices
# end



function global_to_subdomain_bounds(globalarraysize::NTuple{1,T}, topology::CartesianTopology,p) where {T <: Integer}
    ilo, ihi = tile_indices_1d1(globalarraysize[1], topology.global_dims[1], topology.coords[1])
    return (ilo, ihi)
end

function global_to_subdomain_bounds(globalarraysize::NTuple{2,T}, topology::CartesianTopology,p) where {T <: Integer}
    # @show globalarraysize
    # @show topology.global_dims
    # @show topology.coords

    topo_coords = MPI.Cart_coords(topology.comm) |> reverse
    # @show topo_coords, topology.rank

    # @error "adfadsf"
    ilo, ihi = tile_indices_1d1(globalarraysize[1], topology.global_dims[1], topology.coords[1])
    jlo, jhi = tile_indices_1d1(globalarraysize[2], topology.global_dims[2], topology.coords[2])
    # @show (ilo, ihi, jlo, jhi)
    # MPI.Abort(topology.comm, 10)
    return (ilo, ihi, jlo, jhi)
end

function global_to_subdomain_bounds(globalarraysize::NTuple{3,T}, topology::CartesianTopology,p) where {T <: Integer}
    ilo, ihi = tile_indices_1d1(globalarraysize[1], topology.global_dims[1], topology.coords[1])
    jlo, jhi = tile_indices_1d1(globalarraysize[2], topology.global_dims[2], topology.coords[2])
    klo, khi = tile_indices_1d1(globalarraysize[3], topology.global_dims[3], topology.coords[3])
    return (ilo, ihi, jlo, jhi, klo, khi)
end

# """Given tile id in a 1D layout, returns the corresponding tile indices in a 2D layout"""
# function tile_id_to_ij(id::Integer, ntiles::Integer)
#     if id < 1
#         @error("Invalid tile id")
#     end
#     I, J = num_2d_tiles(ntiles)
#     CI = CartesianIndices((1:I, 1:J))
#     ij = Tuple(CI[id])
#     return ij
# end

# """Given tile id in a 1D layout, returns the corresponding tile indices in a 3D layout"""
# function tile_id_to_ijk(id::Integer, ntiles::Integer)
#     if id < 1
#         @error("Invalid tile id")
#     end
#     I, J, K = num_3d_tiles(ntiles)
#     CI = CartesianIndices((1:I, 1:J, 1:K))
#     ijk = Tuple(CI[id])
#     return ijk
# end

# end
