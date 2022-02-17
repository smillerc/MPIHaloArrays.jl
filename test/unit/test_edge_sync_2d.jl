using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This MPIHaloArray edge_sync test is designed with 4 processes only"

function test_edge_sync_2darray_2halo_no_periodic()
    topology = CartesianTopology([4,2], [false, false])

    nhalo = 2
    ni = 8
    nj = 8

    A = MPIHaloArray(zeros(Int, ni,nj), topology, nhalo)
    # Fill each by the current rank
    filldomain!(A, A.topology.rank)
    fillhalo!(A, -1)

    ilo, ihi  = A.local_indices[1].domain
    jlo, jhi  = A.local_indices[2].domain

    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo

    A_ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi]
    A_ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi]
    A_jlo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end]
    A_jhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end]

    sync_edges!(A)

    if rank == 0
        @test all(A_ilo_halo .== -1)
        @test all(A_ihi_halo .== 1)
        @test all(A_jlo_halo .== -1)
        @test all(A_jhi_halo .== 4)
    elseif rank == 1
        @test all(A_ilo_halo .== 0)
        @test all(A_ihi_halo .== 2)
        @test all(A_jlo_halo .== -1)
        @test all(A_jhi_halo .== 5)
    elseif rank == 2
        @test all(A_ilo_halo .== 1)
        @test all(A_ihi_halo .== 3)
        @test all(A_jlo_halo .== -1)
        @test all(A_jhi_halo .== 6)
    elseif rank == 3
        @test all(A_ilo_halo .== 2)
        @test all(A_ihi_halo .== -1)
        @test all(A_jlo_halo .== -1)
        @test all(A_jhi_halo .== 7)
    elseif rank == 4
        @test all(A_ilo_halo .== -1)
        @test all(A_ihi_halo .== 5)
        @test all(A_jlo_halo .== 0)
        @test all(A_jhi_halo .== -1)
    elseif rank == 5
        @test all(A_ilo_halo .== 4)
        @test all(A_ihi_halo .== 6)
        @test all(A_jlo_halo .== 1)
        @test all(A_jhi_halo .== -1)
    elseif rank == 6
        @test all(A_ilo_halo .== 5)
        @test all(A_ihi_halo .== 7)
        @test all(A_jlo_halo .== 2)
        @test all(A_jhi_halo .== -1)
    elseif rank == 7
        @test all(A_ilo_halo .== 6)
        @test all(A_ihi_halo .== -1)
        @test all(A_jlo_halo .== 3)
        @test all(A_jhi_halo .== -1)
    end
end

function test_edge_sync_2darray_2halo_all_periodic()
    topology = CartesianTopology([4,2], [true, true])

    nhalo = 2
    ni = 8
    nj = 8

    A = MPIHaloArray(zeros(Int, ni,nj), topology, nhalo)
    # Fill each by the current rank
    filldomain!(A, A.topology.rank)
    fillhalo!(A, -1)

    ilo, ihi  = A.local_indices[1].domain
    jlo, jhi  = A.local_indices[2].domain

    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo

    A_ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi]
    A_ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi]
    A_jlo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end]
    A_jhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end]

    sync_edges!(A)

    if rank == 0
        @test all(A_ilo_halo .== 3)
        @test all(A_ihi_halo .== 1)
        @test all(A_jlo_halo .== 4)
        @test all(A_jhi_halo .== 4)
    elseif rank == 1
        @test all(A_ilo_halo .== 0)
        @test all(A_ihi_halo .== 2)
        @test all(A_jlo_halo .== 5)
        @test all(A_jhi_halo .== 5)
    elseif rank == 2
        @test all(A_ilo_halo .== 1)
        @test all(A_ihi_halo .== 3)
        @test all(A_jlo_halo .== 6)
        @test all(A_jhi_halo .== 6)
    elseif rank == 3
        @test all(A_ilo_halo .== 2)
        @test all(A_ihi_halo .== 0)
        @test all(A_jlo_halo .== 7)
        @test all(A_jhi_halo .== 7)
    elseif rank == 4
        @test all(A_ilo_halo .== 7)
        @test all(A_ihi_halo .== 5)
        @test all(A_jlo_halo .== 0)
        @test all(A_jhi_halo .== 0)
    elseif rank == 5
        @test all(A_ilo_halo .== 4)
        @test all(A_ihi_halo .== 6)
        @test all(A_jlo_halo .== 1)
        @test all(A_jhi_halo .== 1)
    elseif rank == 6
        @test all(A_ilo_halo .== 5)
        @test all(A_ihi_halo .== 7)
        @test all(A_jlo_halo .== 2)
        @test all(A_jhi_halo .== 2)
    elseif rank == 7
        @test all(A_ilo_halo .== 6)
        @test all(A_ihi_halo .== 4)
        @test all(A_jlo_halo .== 3)
        @test all(A_jhi_halo .== 3)
    end
end

test_edge_sync_2darray_2halo_no_periodic()
test_edge_sync_2darray_2halo_all_periodic()

GC.gc()
MPI.Finalize()
