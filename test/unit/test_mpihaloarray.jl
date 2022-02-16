using MPI
using MPIHaloArrays
# include("../src/MPIHaloArrays.jl")

MPI.Init()

# Create the MPI topology, which here is a Cartesian 4x4 domain (using 16 cores)
topology = CartesianTopology([4,4], [true, true])
rank = MPI.Comm_rank(MPI.COMM_WORLD)

# How many halo cells in each dimension (fixed for all dimensions)
nhalo = 2
ni = 100
nj = 200

# Create an uninitialized array with the given topology
# x = MPIHaloArray{Float64}(ni, nj, topo, nhalo)
x = MPIHaloArray(zeros(ni,nj), topology, nhalo)

# Fill each by the current rank
fill!(x, rank)

GC.gc()
MPI.Finalize()