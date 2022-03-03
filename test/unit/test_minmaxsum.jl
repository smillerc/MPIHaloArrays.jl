using MPI, Test
using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)
const root = 0

@assert nprocs == 8 "This MPIHaloArray edge_sync test is designed with 4 processes only"

function print_array(U,proc)
    if rank == proc
        println("proc: ", proc)
        display(U.data)
        println()
    end
    MPI.Barrier(comm)
end

function test_2d_minmax()
    topology = CartesianTopology(comm, [4,2], [false, false])

    nhalo = 2
    ni = 6
    nj = 5
    data = ones(ni,nj) * rank
    A = MPIHaloArray(data, topology, nhalo)
    fillhalo!(A, 0)

    gmax = globalmax(A)
    gmin = globalmin(A)
    if rank == root
        @test gmax == nprocs - 1
        @test gmin == 0
    end

    gmax = globalmax(A; broadcast=true)
    gmin = globalmin(A; broadcast=true)
    @test gmax == nprocs - 1
    @test gmin == 0
end

function test_2d_sum()
    topology = CartesianTopology(comm, [4,2], [false, false])

    nhalo = 2
    ni = 6
    nj = 5
    data = ones(ni,nj) * rank
    actual_global_sum = sum(sum(ones(ni,nj)*r) for r in 0:nprocs-1)
    A = MPIHaloArray(data, topology, nhalo)
    fillhalo!(A, 0)

    gsum = globalsum(A)
    if rank == root
        @test gsum == actual_global_sum
    end

    gsum = globalsum(A; broadcast=true)
    @test gsum == actual_global_sum
end

test_2d_minmax()
test_2d_sum()