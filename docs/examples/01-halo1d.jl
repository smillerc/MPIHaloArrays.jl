using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This example is designed with 8 processes..."

function print_haloarray(A)
    for p in 0:nprocs-1
        if rank == p
            println("Rank $(p): $(A)")
        end
        MPI.Barrier(comm)
    end
end

topology = CartesianTopology(8, false)

nhalo = 2
ni = 8
data = collect(1:ni) * (rank + 10)
A = MPIHaloArray(data, topology, nhalo)
fillhalo!(A, -1)

if rank == 0 println("Before Sync") end
print_haloarray(A)

updatehalo!(A)

if rank == 0 println("\nAfter Sync") end
print_haloarray(A)

if rank == 0 println("Note that the low boundary on 0 and high boundary on 7 are -1 (non-periodic)") end
GC.gc()
MPI.Finalize()
