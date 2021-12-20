using MPI, Test
using InteractiveUtils

include("../src/topology.jl")

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 16 "Topology neighbor tests are designed with 16 processes only"

function test_4x4_topology()
    P = CartesianTopology([4,4], [true, true])

    # Layout for 4x4 domain; rank and (coords) are shown
    #   |-------|-------|-------|-------|
    #   |       |       |       |       |
    #   |   0   |   1   |   2   |   3   |
    #   | (0,0) | (0,1) | (0,2) | (0,3) |
    #   |-------|-------|-------|-------|
    #   |       |       |       |       |
    #   |   4   |   5   |   6   |   7   |
    #   | (1,0) | (1,1) | (1,2) | (1,3) |
    #   |-------|-------|-------|-------|
    #   |       |       |       |       |
    #   |   8   |   9   |   10  |   11  |
    #   | (2,0) | (2,1) | (2,2) | (2,3) |
    #   |-------|-------|-------|-------|
    #   |       |       |       |       |
    #   |   12  |   13  |   14  |   15  |
    #   | (3,0) | (3,1) | (3,2) | (3,3) |
    #   |-------|-------|-------|-------|

    if P.rank == 0 # test a corner
        @test ilo_neighbor(P) == 3
        @test ihi_neighbor(P) == 1
        @test jlo_neighbor(P) == 4
        @test jhi_neighbor(P) == 12
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1

        @test neighbor(P,1,1,0) == 13
        @test neighbor(P,1,-1,0) == 5  
        @test neighbor(P,-1,1,0) == 15  
        @test neighbor(P,-1,-1,0) == 7  

    elseif P.rank == 5 # test a middle node
        @test ilo_neighbor(P) == 4
        @test ihi_neighbor(P) == 6
        @test jlo_neighbor(P) == 9
        @test jhi_neighbor(P) == 1
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1
        
        @test neighbor(P,1,1,0) == 2  
        @test neighbor(P,1,-1,0) == 10  
        @test neighbor(P,-1,1,0) == 0  
        @test neighbor(P,-1,-1,0) == 8  

    elseif P.rank == 13 # test one on the edge
        @test ilo_neighbor(P) == 12
        @test ihi_neighbor(P) == 14
        @test jlo_neighbor(P) == 1
        @test jhi_neighbor(P) == 9
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1
        
        @test neighbor(P,1,1,0) == 10  
        @test neighbor(P,1,-1,0) == 2  
        @test neighbor(P,-1,1,0) == 8  
        @test neighbor(P,-1,-1,0) == 0  
    end
end

# test_4x4_topology()

GC.gc()
MPI.Finalize()