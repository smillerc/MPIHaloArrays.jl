using MPI, Test
using InteractiveUtils
using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

const ilo = -1
const ihi = +1
const jlo = -1
const jhi = +1
const klo = -1
const khi = +1
const i = 0
const j = 0
const k = 0

@assert nprocs == 16 "Topology neighbor tests are designed with 16 processes only"

function test_4x4_topology_all_periodic()
    P = CartesianTopology([4,4], [true, true])

    # Layout for 4x4 domain; rank and (coords) are shown
    
    
    #   +----- (i)
    #   |
    #   | 
    #   |  (j)

    #      ilo --> ihi
    #   |-------|-------|-------|-------|   jlo
    #   |       |       |       |       |    |
    #   |   0   |   1   |   2   |   3   |    ↓   
    #   | (0,0) | (0,1) | (0,2) | (0,3) |   jhi
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
        @test jlo_neighbor(P) == 12
        @test jhi_neighbor(P) == 4
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1

        @test neighbor(P, ihi, jhi) == 5
        @test neighbor(P, ihi, jlo) == 13  
        @test neighbor(P, ilo, jhi) == 7  
        @test neighbor(P, ilo, jlo) == 15  

    elseif P.rank == 5 # test a middle node
        @test ilo_neighbor(P) == 4
        @test ihi_neighbor(P) == 6
        @test jlo_neighbor(P) == 1
        @test jhi_neighbor(P) == 9
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1
        
        @test neighbor(P, ihi, jhi) == 10
        @test neighbor(P, ihi, jlo) == 2
        @test neighbor(P, ilo, jhi) == 8  
        @test neighbor(P, ilo, jlo) == 0  

    elseif P.rank == 13 # test one on the edge
        @test ilo_neighbor(P) == 12
        @test ihi_neighbor(P) == 14
        @test jlo_neighbor(P) == 9
        @test jhi_neighbor(P) == 1
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1
        
        @test neighbor(P, ihi, jhi) == 2
        @test neighbor(P, ihi, jlo) == 10
        @test neighbor(P, ilo, jhi) == 0
        @test neighbor(P, ilo, jlo) == 8  
    end
end

function test_4x4_topology_no_periodic()
    P = CartesianTopology([4,4], [false, false])

    # Layout for 4x4 domain; rank and (coords) are shown
    
    
    #   +----- (i)
    #   |
    #   | 
    #   |  (j)

    #      ilo --> ihi
    #   |-------|-------|-------|-------|   jlo
    #   |       |       |       |       |    |
    #   |   0   |   1   |   2   |   3   |    ↓   
    #   | (0,0) | (0,1) | (0,2) | (0,3) |   jhi
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
        @test ilo_neighbor(P) == -1
        @test ihi_neighbor(P) == 1
        @test jlo_neighbor(P) == -1
        @test jhi_neighbor(P) == 4
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1

        @test neighbor(P, ihi, jhi) == 5
        @test neighbor(P, ihi, jlo) == -1  
        @test neighbor(P, ilo, jhi) == -1 
        @test neighbor(P, ilo, jlo) == -1  

    elseif P.rank == 5 # test a middle node
        @test ilo_neighbor(P) == 4
        @test ihi_neighbor(P) == 6
        @test jlo_neighbor(P) == 1
        @test jhi_neighbor(P) == 9
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1
        
        @test neighbor(P, ihi, jhi) == 10
        @test neighbor(P, ihi, jlo) == 2
        @test neighbor(P, ilo, jhi) == 8  
        @test neighbor(P, ilo, jlo) == 0  

    elseif P.rank == 13 # test one on the edge
        @test ilo_neighbor(P) == 12
        @test ihi_neighbor(P) == 14
        @test jlo_neighbor(P) == 9
        @test jhi_neighbor(P) == -1
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) == -1
        
        @test neighbor(P, ihi, jhi) == -1
        @test neighbor(P, ihi, jlo) == 10
        @test neighbor(P, ilo, jhi) == -1
        @test neighbor(P, ilo, jlo) == 8  
    end
end

test_4x4_topology_all_periodic()
test_4x4_topology_no_periodic()

GC.gc()
MPI.Finalize()