using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This MPIHaloArray edge_sync test is designed with 4 processes only"

function test_edge_sync_3darray_2halo_no_periodic()
    topology = CartesianTopology([2,2,2], [false,false,false])

    nhalo = 2
    ni = 8
    nj = 8
    nk = 8

    A = MPIHaloArray(zeros(Int, ni, nj, nk), topology, nhalo; do_corners=true)
    # A = MPIHaloArray(zeros(Int, ni), topology, nhalo)
    # Fill each by the current rank
    filldomain!(A, A.topology.rank)
    fillhalo!(A, -1)

    ilo, ihi  = A.local_indices[1].domain
    jlo, jhi  = A.local_indices[2].domain
    klo, khi  = A.local_indices[3].domain

    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    klo_halo_start, klo_halo_end = A.local_indices[3].lo_halo
    khi_halo_start, khi_halo_end = A.local_indices[3].hi_halo

    A_ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, klo:khi]
    A_ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, klo:khi]
    A_jlo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, klo:khi]
    A_jhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, klo:khi]
    A_klo_halo = @view A.data[ilo:ihi, jlo:jhi, klo_halo_start:klo_halo_end]
    A_khi_halo = @view A.data[ilo:ihi, jlo:jhi, khi_halo_start:khi_halo_end]

    A_ilojlo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo:khi]
    A_ilojhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo:khi]
    A_ihijlo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo:khi]
    A_ihijhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo:khi]
    
    A_jloklo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
    A_jlokhi_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
    A_jhiklo_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
    A_jhikhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]
    
    A_iloklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, klo_halo_start:klo_halo_end]
    A_ilokhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, khi_halo_start:khi_halo_end]
    A_ihiklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, klo_halo_start:klo_halo_end]
    A_ihikhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, khi_halo_start:khi_halo_end]

    A_ilojloklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
    A_ilojhiklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
    A_ihijloklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
    A_ihijhiklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
    
    A_ilojlokhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
    A_ilojhikhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]
    A_ihijlokhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
    A_ihijhikhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]

    updatehalo!(A)

    if rank == 0
        @test all(A_ilo_halo .== -1)
        @test all(A_ihi_halo .== 1)
        @test all(A_jlo_halo .== -1)
        @test all(A_jhi_halo .== 2)
        @test all(A_klo_halo .== -1)
        @test all(A_khi_halo .== 4)

        @test all(A_ilojlo_halo .== -1)
        @test all(A_ilojhi_halo .== -1)
        @test all(A_ihijlo_halo .== -1)
        @test all(A_ihijhi_halo .== 3)

        @test all(A_jloklo_halo .== -1)
        @test all(A_jlokhi_halo .== -1)
        @test all(A_jhiklo_halo .== -1)
        @test all(A_jhikhi_halo .== 6)

        @test all(A_iloklo_halo .== -1)
        @test all(A_ilokhi_halo .== -1)
        @test all(A_ihiklo_halo .== -1)
        @test all(A_ihikhi_halo .== 5)

        @test all(A_ilojloklo_halo .== -1)
        @test all(A_ilojhiklo_halo .== -1)
        @test all(A_ihijloklo_halo .== -1)
        @test all(A_ihijhiklo_halo .== -1)
        
        @test all(A_ilojlokhi_halo .== -1)
        @test all(A_ilojhikhi_halo .== -1)
        @test all(A_ihijlokhi_halo .== -1)
        @test all(A_ihijhikhi_halo .== 7)
    end
end

function test_edge_sync_3darray_2halo_all_periodic()
    topology = CartesianTopology([2,2,2], [true, true, true])

    nhalo = 2
    ni = 8
    nj = 8
    nk = 8

    A = MPIHaloArray(zeros(Int, ni, nj, nk), topology, nhalo; do_corners=true)
    # A = MPIHaloArray(zeros(Int, ni), topology, nhalo)
    # Fill each by the current rank
    filldomain!(A, A.topology.rank)
    fillhalo!(A, -1)

    ilo, ihi  = A.local_indices[1].domain
    jlo, jhi  = A.local_indices[2].domain
    klo, khi  = A.local_indices[3].domain

    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    klo_halo_start, klo_halo_end = A.local_indices[3].lo_halo
    khi_halo_start, khi_halo_end = A.local_indices[3].hi_halo

    A_ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, klo:khi]
    A_ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, klo:khi]
    A_jlo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, klo:khi]
    A_jhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, klo:khi]
    A_klo_halo = @view A.data[ilo:ihi, jlo:jhi, klo_halo_start:klo_halo_end]
    A_khi_halo = @view A.data[ilo:ihi, jlo:jhi, khi_halo_start:khi_halo_end]

    A_ilojlo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo:khi]
    A_ilojhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo:khi]
    A_ihijlo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo:khi]
    A_ihijhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo:khi]
    
    A_jloklo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
    A_jlokhi_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
    A_jhiklo_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
    A_jhikhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]
    
    A_iloklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, klo_halo_start:klo_halo_end]
    A_ilokhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, khi_halo_start:khi_halo_end]
    A_ihiklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, klo_halo_start:klo_halo_end]
    A_ihikhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, khi_halo_start:khi_halo_end]

    A_ilojloklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
    A_ilojhiklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
    A_ihijloklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
    A_ihijhiklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
    
    A_ilojlokhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
    A_ilojhikhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]
    A_ihijlokhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
    A_ihijhikhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]

    updatehalo!(A)

    if rank == 0
        @test all(A_ilo_halo .== 1)
        @test all(A_ihi_halo .== 1)
        @test all(A_jlo_halo .== 2)
        @test all(A_jhi_halo .== 2)
        @test all(A_klo_halo .== 4)
        @test all(A_khi_halo .== 4)

        @test all(A_ilojlo_halo .== 3)
        @test all(A_ilojhi_halo .== 3)
        @test all(A_ihijlo_halo .== 3)
        @test all(A_ihijhi_halo .== 3)

        @test all(A_jloklo_halo .== 6)
        @test all(A_jlokhi_halo .== 6)
        @test all(A_jhiklo_halo .== 6)
        @test all(A_jhikhi_halo .== 6)

        @test all(A_iloklo_halo .== 5)
        @test all(A_ilokhi_halo .== 5)
        @test all(A_ihiklo_halo .== 5)
        @test all(A_ihikhi_halo .== 5)

        @test all(A_ilojloklo_halo .== 7)
        @test all(A_ilojhiklo_halo .== 7)
        @test all(A_ihijloklo_halo .== 7)
        @test all(A_ihijhiklo_halo .== 7)

        @test all(A_ilojlokhi_halo .== 7)
        @test all(A_ilojhikhi_halo .== 7)
        @test all(A_ihijlokhi_halo .== 7)
        @test all(A_ihijhikhi_halo .== 7)
    end
end


test_edge_sync_3darray_2halo_no_periodic()
test_edge_sync_3darray_2halo_all_periodic()

GC.gc()
MPI.Finalize()
