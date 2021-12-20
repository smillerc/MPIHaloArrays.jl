using MPI, Test
using InteractiveUtils

include("../src/topology.jl")

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 16 "Topology neighbor tests are designed with 16 processes only"

function test_2x2x4_topology()
    P = CartesianTopology([2,2,4], [true, true, true])

    # Layout for 4x4 domain; rank and (coords) are shown
    # ordering 
    # i: left to right
    # j: top to bottom
    # k: back to front; k = 0 is the back-most slab

    #   |---------|---------|  ..  |---------|---------|  .. |---------|---------|  .. |---------|---------| 
    #   |         |         |  ..  |         |         |  .. |         |         |  .. |         |         | 
    #   |    0    |    1    |  ..  |    4    |    5    |  .. |    8    |    9    |  .. |    12   |    13   | 
    #   | (0,0,0) | (0,1,0) |  ..  | (0,0,1) | (0,1,1) |  .. | (0,0,1) | (0,1,2) |  .. | (0,0,3) | (0,1,3) | 
    #   |---------|---------|  ..  |---------|---------|  .. |---------|---------|  .. |---------|---------| 
    #   |         |         |  ..  |         |         |  .. |         |         |  .. |         |         | 
    #   |    2    |    3    |  ..  |    6    |    7    |  .. |    10   |    11   |  .. |    14   |    15   | 
    #   | (1,0,0) | (1,1,0) |  ..  | (1,0,1) | (1,1,1) |  .. | (1,0,1) | (1,1,2) |  .. | (1,0,3) | (1,1,3) | 
    #   |---------|---------|  ..  |---------|---------|  .. |---------|---------|  .. |---------|---------| 
    #      "back most", k=0                                                           "front most" slab, k=3

    # test a corner located at mpi coords (0,0,0), 
    # but at (imin, jmax, kmax) in array coords
    if P.rank == 0 
        @test ilo_neighbor(P) == 1
        @test ihi_neighbor(P) == 1
        @test jlo_neighbor(P) == 2
        @test jhi_neighbor(P) == 2
        @test klo_neighbor(P) == 4
        @test khi_neighbor(P) == 12

        @test neighbor(P, 1, 1,0) == 3
        @test neighbor(P, 1,-1,0) == 3
        @test neighbor(P,-1,-1,0) == 3
        @test neighbor(P,-1, 1,0) == 3

        @test neighbor(P,0,0,1) == 12
        @test neighbor(P,1,0,1) == 13
        # @test neighbor(P,1,-1,0) == 5  
        # @test neighbor(P,-1,1,0) == 15  
        # @test neighbor(P,-1,-1,0) == 7  

    # elseif P.rank == 5 # test a middle node
    #     @test ilo_neighbor(P) == 4
    #     @test ihi_neighbor(P) == 6
    #     @test jlo_neighbor(P) == 9
    #     @test jhi_neighbor(P) == 1
    #     @test klo_neighbor(P) == -1
    #     @test khi_neighbor(P) == -1
        
    #     @test neighbor(P,1,1,0) == 2  
    #     @test neighbor(P,1,-1,0) == 10  
    #     @test neighbor(P,-1,1,0) == 0  
    #     @test neighbor(P,-1,-1,0) == 8  

    # elseif P.rank == 13 # test one on the edge
    #     @test ilo_neighbor(P) == 12
    #     @test ihi_neighbor(P) == 14
    #     @test jlo_neighbor(P) == 1
    #     @test jhi_neighbor(P) == 9
    #     @test klo_neighbor(P) == -1
    #     @test khi_neighbor(P) == -1
        
    #     @test neighbor(P,1,1,0) == 10  
    #     @test neighbor(P,1,-1,0) == 2  
    #     @test neighbor(P,-1,1,0) == 8  
    #     @test neighbor(P,-1,-1,0) == 0  
    end
end



test_2x2x4_topology()

GC.gc()
MPI.Finalize()