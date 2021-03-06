using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This example is designed with 8 processes..."

topology = CartesianTopology(comm, [2,2,2], [false,false,false])

nhalo = 2
ni = 8
nj = 8
nk = 8

A = MPIHaloArray(zeros(Int, ni, nj, nk), topology, nhalo; do_corners=true)
filldomain!(A, rank)
fillhalo!(A, -1)

updatehalo!(A)

GC.gc()
MPI.Finalize()
