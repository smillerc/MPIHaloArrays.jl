using MPI, Test
using OffsetArrays
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 16 "MPIHaloArray tests are designed with 16 processes only"

function test_2d_indices()
    topology = CartesianTopology([4,4], [false, false])
    # topology = CartesianTopology([1], [false])

    nhalo = 2
    ni = 8
    nj = 6
    A = MPIHaloArray(zeros(Int, ni,nj), topology, nhalo)

    ilo_halo_start, ilo_halo_end, ilo_dom_start, ilo_dom_end = lo_indices(A, 1)
    jlo_halo_start, jlo_halo_end, jlo_dom_start, jlo_dom_end = lo_indices(A, 2)
    ihi_dom_start, ihi_dom_end, ihi_halo_start, ihi_halo_end = hi_indices(A, 1)
    jhi_dom_start, jhi_dom_end, jhi_halo_start, jhi_halo_end = hi_indices(A, 2)

    @test ilo_halo_start == 1
    @test ilo_halo_end == 2
    @test ilo_dom_start == 3
    @test ilo_dom_end == 4

    @test jlo_halo_start == 1
    @test jlo_halo_end == 2
    @test jlo_dom_start == 3
    @test jlo_dom_end == 4

    # 1 2 | 3 4 5 6 | 7 8
    @test ihi_dom_start == 5
    @test ihi_dom_end == 6
    @test ihi_halo_start == 7
    @test ihi_halo_end == 8

    # 1 2 | 3 4 |5 6
    @test jhi_dom_start == 3
    @test jhi_dom_end == 4
    @test jhi_halo_start == 5
    @test jhi_halo_end == 6
end

function test_2d_indices_offset()
    topology = CartesianTopology([4,4], [false, false])

    nhalo = 2
    offsetdata = OffsetArray(rand(Int,6,8), -2:3, -4:3)
    A = MPIHaloArray(offsetdata, topology, nhalo)

    ilo_halo_start, ilo_halo_end, ilo_dom_start, ilo_dom_end = lo_indices(A, 1)
    jlo_halo_start, jlo_halo_end, jlo_dom_start, jlo_dom_end = lo_indices(A, 2)
    ihi_dom_start, ihi_dom_end, ihi_halo_start, ihi_halo_end = hi_indices(A, 1)
    jhi_dom_start, jhi_dom_end, jhi_halo_start, jhi_halo_end = hi_indices(A, 2)

    # -2 -1 | 0 1 | 2 3
    @test ilo_halo_start == -2
    @test ilo_halo_end == -1
    @test ilo_dom_start == 0
    @test ilo_dom_end == 1
    @test ihi_dom_start == 0
    @test ihi_dom_end == 1
    @test ihi_halo_start == 2
    @test ihi_halo_end == 3

    # -4 -3 | -2 -1 0 1 | 2 3
    @test jlo_halo_start == -4
    @test jlo_halo_end == -3
    @test jlo_dom_start == -2
    @test jlo_dom_end == -1
    @test jhi_dom_start == 0
    @test jhi_dom_end == 1
    @test jhi_halo_start == 2
    @test jhi_halo_end == 3
end

test_2d_indices()
test_2d_indices_offset()

GC.gc()
MPI.Finalize()