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

function test_2x2x4_topology_all_periodic()
    P = CartesianTopology(comm, [2,2,4], [true, true, true])

    # klo
    # 15   14  15   jhi
    # 13   12  13   j
    # 15   14  15   jlo
    # ilo  i   ihi
    
    # k
    #  3   2   3   jhi
    #  1   0   1   j
    #  3   2   3   jlo
    # ilo  i  ihi

    # khi
    #  7   6   7   jhi
    #  5   4   5   j
    #  7   6   7   jlo
    # ilo  i  ihi

    if P.rank == 0 
        # klo
        @test neighbor(P, ilo, jlo, klo) == 15
        @test neighbor(P, i  , jlo, klo) == 14
        @test neighbor(P, ihi, jlo, klo) == 15
        
        @test neighbor(P, ilo, j  , klo) == 13
        @test neighbor(P, ihi, j  , klo) == 13

        @test neighbor(P, ilo, jhi, klo) == 15
        @test neighbor(P, i  , jhi, klo) == 14
        @test neighbor(P, ihi, jhi, klo) == 15

        # k
        @test ilo_neighbor(P) == 1
        @test ihi_neighbor(P) == 1
        @test jlo_neighbor(P) == 2
        @test jhi_neighbor(P) == 2
        @test klo_neighbor(P) == 12
        @test khi_neighbor(P) == 4

        @test neighbor(P, ihi, jhi, k) == 3
        @test neighbor(P, ihi, jlo, k) == 3
        @test neighbor(P, ilo, jlo, k) == 3
        @test neighbor(P, ilo, jhi, k) == 3

        # khi
        @test neighbor(P, ilo, jlo, khi) == 7
        @test neighbor(P, i  , jlo, khi) == 6
        @test neighbor(P, ihi, jlo, khi) == 7
        @test neighbor(P, ilo, j  , khi) == 5
        @test neighbor(P, ihi, j  , khi) == 5
        @test neighbor(P, ilo, jhi, khi) == 7
        @test neighbor(P, i  , jhi, khi) == 6
        @test neighbor(P, ihi, jhi, khi) == 7
    end
end

function test_2x2x4_topology_no_periodic()
    P = CartesianTopology(comm, [2,2,4], [false, false, false])


    #      / (k)
    #    /
    #   +----- (i)
    #   |
    #   | (j)
    #
    
    # klo
    #  -1   -1  -1  jhi
    #  -1   -1  -1  j
    #  -1   -1  -1   jlo
    # ilo  i   ihi

    # k
    #  -1  -1  -1   jhi
    #  -1   0   1   j
    #  -1   2   3   jlo
    # ilo  i  ihi

    # khi
    #  -1  -1  -1   jhi
    #  -1   4   5   j
    #  -1   6   7   jlo
    # ilo  i  ihi

    if P.rank == 0 
        # klo
        @test neighbor(P, ilo, jlo, klo) == -1
        @test neighbor(P, i  , jlo, klo) == -1
        @test neighbor(P, ihi, jlo, klo) == -1
        @test neighbor(P, ilo, j  , klo) == -1
        @test neighbor(P, ihi, j  , klo) == -1
        @test neighbor(P, ilo, jhi, klo) == -1
        @test neighbor(P, i  , jhi, klo) == -1
        @test neighbor(P, ihi, jhi, klo) == -1

        # k
        @test ilo_neighbor(P) == -1
        @test ihi_neighbor(P) ==  1
        @test jlo_neighbor(P) == -1
        @test jhi_neighbor(P) ==  2
        @test klo_neighbor(P) == -1
        @test khi_neighbor(P) ==  4

        @test neighbor(P, ihi, jhi, k) == 3
        @test neighbor(P, ihi, jlo, k) == -1
        @test neighbor(P, ilo, jlo, k) == -1
        @test neighbor(P, ilo, jhi, k) == -1

        # khi
        @test neighbor(P, ilo, jlo, khi) == -1
        @test neighbor(P, i  , jlo, khi) == -1
        @test neighbor(P, ihi, jlo, khi) == -1
        @test neighbor(P, ilo, j  , khi) == -1
        @test neighbor(P, ihi, j  , khi) == 5
        @test neighbor(P, ilo, jhi, khi) == -1
        @test neighbor(P, i  , jhi, khi) == 6
        @test neighbor(P, ihi, jhi, khi) == 7
    end
end

test_2x2x4_topology_all_periodic()
test_2x2x4_topology_no_periodic()

GC.gc()
MPI.Finalize()