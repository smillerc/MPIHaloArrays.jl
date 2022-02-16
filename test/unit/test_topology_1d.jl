using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

const ilo = -1
const ihi = +1
const i = 0

@assert nprocs == 16 "Topology neighbor tests are designed with 16 processes only"

function test_16x1_topology_all_periodic()
    P = CartesianTopology([16], [true])

    #  0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15
    if P.rank == 0 # test at edge
        @test ilo_neighbor(P) == 15
        @test ihi_neighbor(P) == 1

    elseif P.rank == 5 # test a middle node
        @test ilo_neighbor(P) == 4
        @test ihi_neighbor(P) == 6
    end
end

function test_16x1_topology_no_periodic()
    P = CartesianTopology([16], [false])

    #  0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15
    if P.rank == 0 # test at edge
        @test ilo_neighbor(P) == -1
        @test ihi_neighbor(P) == 1

    elseif P.rank == 5 # test a middle node
        @test ilo_neighbor(P) == 4
        @test ihi_neighbor(P) == 6
    end
end

test_16x1_topology_all_periodic()
test_16x1_topology_no_periodic()

GC.gc()
MPI.Finalize()