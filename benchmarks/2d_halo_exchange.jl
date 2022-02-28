using MPI, MPIHaloArrays
using BenchmarkTools

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This example is designed with 8 processes..."

function benchmark_halo_exchange(A)
    updatehalo!(A)
end

topology = CartesianTopology([4,2], [true, true])
nhalo = 2
ni = 512
nj = 512 
data = rand(Float64, ni, nj)
A = MPIHaloArray(data, topology, nhalo)

@time benchmark_halo_exchange(A)

GC.gc()
MPI.Finalize()
