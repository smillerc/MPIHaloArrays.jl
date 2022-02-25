using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This example is designed with 8 processes..."

topology = CartesianTopology([2,2,2], [false, false, true])

# nhalo = 2
# ni = 6
# nj = 5

# data = rank_data(ni, nj, rank)
# A = MPIHaloArray(data, topology, nhalo)
# fillhalo!(A, -1)

# updatehalo!(A)

GC.gc()
MPI.Finalize()
