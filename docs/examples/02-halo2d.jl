using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This example is designed with 8 processes..."

topology = CartesianTopology([4,2], [false, false])

nhalo = 2
ni = 6
nj = 5

A = MPIHaloArray(rand(ni, nj), topology, nhalo)
fillhalo!(A, -1)

updatehalo!(A)

GC.gc()
MPI.Finalize()
