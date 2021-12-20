using MPI
using OffsetArrays

"""An abstract ParallelTopology type that is extended by either a CartesianTopology or GraphTopology (future)"""
abstract type ParallelTopology end

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
 - `neighbors`: OffsetArray{Int}; Neighbor ranks (including corners), indexed as `[[left, center, right], i, j, k]`
"""
struct CartesianTopology <: ParallelTopology
    comm::MPI.Comm
    nprocs::Int
    rank::Int
    coords::Vector{Int}      # [i, j, k]; coordinates in the
    global_dims::Vector{Int} # [i, j, k]; number of domains in each direction
    isperiodic::Vector{Bool} # [i, j, k]; is this dimension periodic?
    neighbors::OffsetArray{Int, 3}   # [[left, center, right], i, j, k]; defaults to -1 if no neighbor
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
#            +-------+-------+-------+
#          /       /       /       / |
#        +-------+-------+-------+   |
#      /       /       /       / |   |
#    /       /       /       /   |   |
#   +-------+-------+-------+    |   +
#   |       |       |       |    | / |
#   | (1,1) | (0,1) | (1,1) |    +   |
#   |       |       |       |  / |   |
#   +-------+-------+-------+/   |   +
#   |       |       |       |    | / |
#   | (-1,0)| (0,0) | (1,2) |    +   |
#   |       |       |       |  / |   |
#   +-------+-------+-------+/   |   +
#   |       |       |       |    | /
#   |(-1,-1)| (0,-1)| (2,2) |    +
#   |       |       |       |  /
#   +-------+-------+-------+/



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
    coords = MPI.Cart_coords(comm_cart)
    neighbors = OffsetArray(-ones(Int8,3,3,3), -1:1, -1:1, -1:1)

    # Using these constants makes it less bug-prone when referencing neighbor indices
    # due to the way MPI does their indexing
    CENTER = 0
    LEFT = -1
    RIGHT = +1
    BOTTOM = +1
    TOP = -1
    FRONT = +1
    BACK = -1

    # MPI convention is (k, j, i), or (z, y, x) which is annoying
    if length(dims) == 1
        left, right = MPI.Cart_shift(comm_cart, 0, 1)
    elseif length(dims) == 2
        top, bottom = MPI.Cart_shift(comm_cart, 0, 1)
        left, right = MPI.Cart_shift(comm_cart, 1, 1)

        topright    = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT, TOP)
        bottomright = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT, BOTTOM)
        topleft     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT, TOP)
        bottomleft  = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT, BOTTOM)

        neighbors[:,  1, 0] = [topleft   , top   , topright]
        neighbors[:,  0, 0] = [left      , rank  , right]
        neighbors[:, -1, 0] = [bottomleft, bottom, bottomright]

    elseif length(dims) == 3
        front, back = MPI.Cart_shift(comm_cart, 0, 1)
        top, bottom = MPI.Cart_shift(comm_cart, 1, 1)
        left, right = MPI.Cart_shift(comm_cart, 2, 1)

        topright    = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT, TOP,    CENTER)
        topleft     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT , TOP,    CENTER)
        bottomright = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT, BOTTOM, CENTER)
        bottomleft  = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT , BOTTOM, CENTER)

        fronttop         = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, CENTER, TOP,    FRONT)
        frontbottom      = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, CENTER, BOTTOM, FRONT)
        frontleft        = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT,  CENTER, FRONT)
        frontright       = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT, CENTER, FRONT)
        fronttopright    = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT,  TOP,    FRONT)
        frontbottomright = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT,  BOTTOM, FRONT)
        fronttopleft     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT,   TOP,    FRONT)
        frontbottomleft  = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT,   BOTTOM, FRONT)

        backtop         = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, CENTER, TOP,    BACK)
        backbottom      = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, CENTER, BOTTOM, BACK)
        backleft        = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT, CENTER, BACK)
        backright       = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity,  RIGHT, CENTER, BACK)
        backtopright    = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT,  TOP,    BACK)
        backbottomright = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, RIGHT,  BOTTOM, BACK)
        backtopleft     = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT,   TOP,    BACK)
        backbottomleft  = offset_coord_to_rank(comm_cart, mpi_dims, mpi_periodicity, LEFT,   BOTTOM, BACK)

        # Center
        neighbors[:,  1, 0] = [topleft   , top   , topright]
        neighbors[:,  0, 0] = [left      , rank  , right]
        neighbors[:, -1, 0] = [bottomleft, bottom, bottomright]

        # Front
        neighbors[:,  1, 1] = [backtopleft   , backtop   , backtopright]
        neighbors[:,  0, 1] = [backleft      , back      , backright]
        neighbors[:, -1, 1] = [backbottomleft, backbottom, backbottomright]

        # Back
        neighbors[:,  1, -1] = [fronttopleft   , fronttop   , fronttopright]
        neighbors[:,  0, -1] = [frontleft      , front      , frontright]
        neighbors[:, -1, -1] = [frontbottomleft, frontbottom, frontbottomright]

    end

    CartesianTopology(comm_cart, nprocs, rank, coords, dims, periodicity, neighbors)
end

"""Helper function to find rank based on 3D offsets"""
function offset_coord_to_rank(comm, dims, periods, i_offset::Int, j_offset::Int, k_offset::Int)
    coords = MPI.Cart_coords(comm) .+ (k_offset, j_offset, i_offset)
    coord_to_rank(comm, dims, periods, coords)
end

"""Helper function to find rank based on 2D offsets"""
function offset_coord_to_rank(comm, dims, periods, i_offset::Int, j_offset::Int)
    coords = MPI.Cart_coords(comm) .+ (j_offset, i_offset)
    coord_to_rank(comm, dims, periods, coords)
end

"""Helper function to find rank based on coordinates"""
function coord_to_rank(comm, dims, periods, coords)

    isvalid = true

    for i in 1:length(dims)
        if (!periods[i] && (coords[i] >= dims[i] || coords[i] < 0))
            isvalid = false
            break
        end
    end

    if isvalid
        rank = MPI.Cart_rank(comm, coords)
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
klo_neighbor(p::CartesianTopology) = p.neighbors[ 0, 0, 1]

"""Neighbor rank in the k+1 direction"""
khi_neighbor(p::CartesianTopology) = p.neighbors[ 0, 0,-1]

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

# Find the upper right corner neighbor (ihi and jhi side)
ihijhi_corner = neighbor(P,+1,+1,0)
```

"""
function neighbor(p::CartesianTopology, i_offset::Int, j_offset::Int, k_offset::Int)
    p.neighbors[i_offset, j_offset, k_offset]
end