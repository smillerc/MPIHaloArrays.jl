using MPI
using OffsetArrays

export ParallelTopology, CartesianTopology
export neighbor, ilo_neighbor, ihi_neighbor, jlo_neighbor, jhi_neighbor, klo_neighbor, khi_neighbor

include("topology.jl")
# const NNEIGHBORS_PER_DIM = 2     # Number of neighbors per dimension (left neighbor + right neighbor).
# const NDIMS_MPI = 3              # Internally, we set the number of dimensions always to 3 for calls to MPI. This ensures a fixed size for MPI coords, neigbors, etc and in general a simple, easy to read code.

# include("partitioning.jl")

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
    data::Array{T,N}
    nhalo::Int
    rank::Int
    topology::CartesianTopology
    # comm::MPI.Comm
    # window::MPI.Win
    # neighbor_ranks::Vector{Int}
    # coords::Vector{Int}
    # MPIHaloArray{T}(sizes::Vararg{<:Integer,N}, nhalo) where {T,N} = MPIHaloArray(Array{T, 2}(undef, sizes...), nhalo, 0...)
end

function MPIHaloArray(A::Array{T,N}, topo::CartesianTopology, nhalo::Int) where {T,N}
    if topo.rank == 0
        @show topo
    end
    MPIHaloArray(A, nhalo, topo.rank, topo)
end

# Required interface overloads to be an AbstractArray
Base.size(A::MPIHaloArray) = size(A.data)
Base.getindex(A::MPIHaloArray{T,N}, i::Int) where {T,N} = getindex(A.data, i)
Base.getindex(A::MPIHaloArray{T,N}, I::Vararg{Int, N}) where {T,N} = getindex(A.data, I...)
Base.setindex!(A::MPIHaloArray{T,N}, v, i::Int) where {T,N} = setindex!(A.data, v, i)
Base.setindex!(A::MPIHaloArray{T,N}, v, I::Vararg{Int, N}) where {T,N} = setindex!(A.data, v..., I...)



# include("utils/indexing.jl")

# include("sync_edges.jl")
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

# """
#     applylocalfunc(f, A::MPIHaloArray)
# Execute the function f on the part of a owned by the current rank. The local part may be modified by f.
# """
# function applylocalfunc!(f, A::MPIHaloArray)
#     MPI.Win_lock(MPI.LOCK_EXCLUSIVE, A.rank, 0, A.window)
#     result = f(A.data)
#     MPI.Win_unlock(A.rank, A.window)
#     return result
# end


# function init_global_grid(nx::Integer, ny::Integer, nz::Integer; 
#                           dimx::Integer=0, dimy::Integer=0, dimz::Integer=0, nhalo=2,
#                           periodx::Integer=0, periody::Integer=0, periodz::Integer=0, 
#                           disp::Integer=1, reorder::Integer=1, comm::MPI.Comm=MPI.COMM_WORLD, init_MPI::Bool=true, quiet::Bool=false)

#     nxyz          = [nx, ny, nz];
#     dims          = [dimx, dimy, dimz];
#     periods       = [periodx, periody, periodz];

#     dims[(nxyz .== 1) .& (dims .== 0)] .= 1;   # Setting any of nxyz to 1, means that the corresponding dimension must also be 1 in the global grid. Thus, the corresponding dims entry must be 1.
    
#     if (init_MPI)  # NOTE: init MPI only, once the input arguments have been checked.
#         if (MPI.Initialized()) error("MPI is already initialized. Set the argument 'init_MPI=false'."); end
#         MPI.Init();
#     else
#         if (!MPI.Initialized()) error("MPI has not been initialized beforehand. Remove the argument 'init_MPI=false'."); end  # Ensure that MPI is always initialized after init_global_grid().
#     end

#     nprocs    = MPI.Comm_size(comm);
#     MPI.Dims_create!(nprocs, dims);
#     comm_cart = MPI.Cart_create(comm, dims, periods, reorder);
#     me        = MPI.Comm_rank(comm_cart);
#     coords    = MPI.Cart_coords(comm_cart);
#     neighbors = fill(MPI.MPI_PROC_NULL, NNEIGHBORS_PER_DIM, NDIMS_MPI);
#     for i = 1:NDIMS_MPI
#         neighbors[:,i] .= MPI.Cart_shift(comm_cart, i - 1, disp);
#     end
#     return me, neighbors, dims, nprocs, coords, comm_cart # The typical use case requires only these variables; the remaining can be obtained calling get_global_grid() if needed.
# end

# # function print_arr(U)
# #     for proc in 0:nprocs
# #         if me == proc
# #             println()
# #             println("proc: ", proc)
# #             for j in size(U, 2):-1:1
# #                 println("j ", j, ":\t", U[:,j])
# #             end
# #         end
# #         MPI.Barrier(comm_cart)
# #     end
# # end

# # ni = 8
# # nj = 10
# # nhalo = 2
# # me, neighbors, dims, nprocs, coords, comm_cart = init_global_grid(ni, nj, 1)

# # ilo_neighbor, jlo_neighbor, klo_neighbor = neighbors[1,:]
# # ihi_neighbor, jhi_neighbor, khi_neighbor = neighbors[2,:]

# # U = zeros(ni, nj)
# # U .= me

# # print_arr(U)

# # println("After")
# # print_arr(U)

# # MPI.Finalize()
