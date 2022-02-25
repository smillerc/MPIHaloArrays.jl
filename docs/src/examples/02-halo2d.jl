# 2D Halo Example

```julia
# examples/02-halo2d.jl
using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This example is designed with 8 processes..."

topology = CartesianTopology([4,2], [false, false])

# nhalo = 2
# ni = 6
# nj = 5

# data = rank_data(ni, nj, rank)
# A = MPIHaloArray(data, topology, nhalo)
# fillhalo!(A, -1)

# updatehalo!(A)

GC.gc()
MPI.Finalize()
```

```
> mpiexecjl -n 8 julia examples/02-halo2d.jl
```
