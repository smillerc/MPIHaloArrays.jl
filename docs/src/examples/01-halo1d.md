# 1D Halo Example

```julia
# examples/01-halo1d.jl
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

if rank == 0 println("Note that the low boundary on 0 and high boundary on 7 are -1 since the domain is not periodic") end
GC.gc()
MPI.Finalize()
```

```
> mpiexecjl -n 8 julia examples/01-halo1d.jl
Before Sync
Rank 0: [-1, -1, 10, 20, 30, 40, 50, 60, 70, 80, -1, -1]
Rank 1: [-1, -1, 11, 22, 33, 44, 55, 66, 77, 88, -1, -1]
Rank 2: [-1, -1, 12, 24, 36, 48, 60, 72, 84, 96, -1, -1]
Rank 3: [-1, -1, 13, 26, 39, 52, 65, 78, 91, 104, -1, -1]
Rank 4: [-1, -1, 14, 28, 42, 56, 70, 84, 98, 112, -1, -1]
Rank 5: [-1, -1, 15, 30, 45, 60, 75, 90, 105, 120, -1, -1]
Rank 6: [-1, -1, 16, 32, 48, 64, 80, 96, 112, 128, -1, -1]
Rank 7: [-1, -1, 17, 34, 51, 68, 85, 102, 119, 136, -1, -1]

After Sync
Rank 0: [-1, -1, 10, 20, 30, 40, 50, 60, 70, 80, 11, 22]
Rank 1: [70, 80, 11, 22, 33, 44, 55, 66, 77, 88, 12, 24]
Rank 2: [77, 88, 12, 24, 36, 48, 60, 72, 84, 96, 13, 26]
Rank 3: [84, 96, 13, 26, 39, 52, 65, 78, 91, 104, 14, 28]
Rank 4: [91, 104, 14, 28, 42, 56, 70, 84, 98, 112, 15, 30]
Rank 5: [98, 112, 15, 30, 45, 60, 75, 90, 105, 120, 16, 32]
Rank 6: [105, 120, 16, 32, 48, 64, 80, 96, 112, 128, 17, 34]
Rank 7: [112, 128, 17, 34, 51, 68, 85, 102, 119, 136, -1, -1]
Note that the low boundary on 0 and high boundary on 7 are -1 since the domain is not periodic
```
