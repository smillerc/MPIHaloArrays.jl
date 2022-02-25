# module ParallelTopologies

using MPI
using OffsetArrays

# export ParallelTopology, CartesianTopology
# export ilo_neighbor, ihi_neighbor, jlo_neighbor, jhi_neighbor, klo_neighbor, khi_neighbor, neighbor, neighbors

"""An abstract ParallelTopology type that is extended by either a CartesianTopology or GraphTopology (future)"""
abstract type ParallelTopology end

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
struct CartesianTopology <: ParallelTopology
    comm::MPI.Comm
    nprocs::Int
    rank::Int
    coords::Vector{Int}      # [i, j, k]; coordinates in the
    global_dims::Vector{Int} # [i, j, k]; number of domains in each direction
    isperiodic::Vector{Bool} # [i, j, k]; is this dimension periodic?
    neighbors::OffsetArray{Int, 3}   # [[ilo, center, ihi], i, j, k]; defaults to -1 if no neighbor
end

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
#   | (1,1) | (0,1) | (1,1) |   / |     +
#   |       |       |       | /   |   / | 
#   +-------+-------+-------+     | /   |
#   |       |       |       |     +     |
#   | (-1,0)| (0,0) | (1,2) |   / |     +
#   |       |       |       | /   |   / 
#   +-------+-------+-------+     | /   
#   |       |       |       |     + 
#   |(-1,-1)| (0,-1)| (2,2) |   /
#   |       |       |       | /
#   +-------+-------+-------+



"""
Create a CartesianTopology type that holds neighbor information, current rank, etc.

# Arguments
 - `dims`: Dimensions of the domain in each direction, e.g. [4,3] means a total of 12 procs, with 4 in x and 3 in y
 - `periodicity`: Vector of bools to set if the domain is periodic along a specific dimension

# Example
```julia

# Create a topology of 4x4 with periodic boundaries in both directions
P = CartesianTopology([4,4], [true, true])
```
"""
function CartesianTopology(dims::Vector{Int}, periodicity::Vector{Bool}; canreorder = false)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    # use reverse(dims) b/c MPI convention in (k, j, i), but the expected order is (i, j, k)
    mpi_dims = reverse(dims)
    mpi_periodicity = reverse(periodicity)

    @assert length(dims) == length(periodicity) "You must specify periodicity (true/false) for each dimension"
    @assert prod(dims) == nprocs "Too many dimensions for the given number of processes"

    comm_cart = MPI.Cart_create(comm, mpi_dims, mpi_periodicity .|> Int, canreorder)
    # coords = MPI.Cart_coords(comm_cart)
    coords = MPI.Cart_coords(comm_cart) |> reverse
    neighbors = OffsetArray(-ones(Int8,3,3,3), -1:1, -1:1, -1:1)

    # MPI convention is (k, j, i), or (z, y, x) which is annoying
    if length(dims) == 1
        ilo, ihi = MPI.Cart_shift(comm_cart, 0, 1)
        
        neighbors[:,  0, 0] = [ilo, rank, ihi]

    elseif length(dims) == 2
        jlo, jhi = MPI.Cart_shift(comm_cart, 0, 1)
        ilo, ihi = MPI.Cart_shift(comm_cart, 1, 1)

        jhi_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JHI)
        jlo_ihi = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, IHI, JLO)
        jhi_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JHI)
        jlo_ilo = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, ILO, JLO)

        neighbors[:,  1, 0] = [jhi_ilo, jhi , jhi_ihi]
        neighbors[:,  0, 0] = [ilo    , rank, ihi]
        neighbors[:, -1, 0] = [jlo_ilo, jlo , jlo_ihi]

    elseif length(dims) == 3
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

    CartesianTopology(comm_cart, nprocs, rank, coords, dims, periodicity, neighbors)
end

function CartesianTopology(dims::Int, periodicity::Bool; canreorder = false)
    CartesianTopology([dims], [periodicity]; canreorder = canreorder)
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
        rank = -1
    end
    rank
end

function Base.show(io::IO, p::ParallelTopology)
    println(io, "ParallelTopology")
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

Find the neighbor rank based on the offesets in (i,j,k).
This follows the traditional array index convention rather than MPI's version, so an `i_offset=1` will shift up in the array indexing.

# Arguments
 - `p` : CartesianTopology type
 - `i_offset`: Offset in the `i` direction
 - `j_offset`: Offset in the `j` direction
 - `k_offset`: Offset in the `k` direction

# Example:
```julia
# Makes a 4x4 domain with periodic boundaries in both dimensions
P = CartesianTopology([4,4], [true, true])

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

# end