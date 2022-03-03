using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This example is designed with 8 processes..."


function print_array(U, proc)
    if rank == proc
        println("rank: ", proc)
        display(U.data)
        println()
    end
    MPI.Barrier(comm)
end

topology = CartesianTopology(comm, [4,2], [false, false])

nhalo = 2
ni = 6
nj = 5

A = MPIHaloArray(zeros(Int, ni, nj), topology, nhalo)
fillhalo!(A, -1)
filldomain!(A, rank)

if rank == 0 println("Before sync") end
for p in 0:nprocs-1
    print_array(A, p)
end

updatehalo!(A)

if rank == 0 println("After sync") end
for p in 0:nprocs-1
    print_array(A, p)
end

GC.gc()
MPI.Finalize()
