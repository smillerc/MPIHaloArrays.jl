using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This MPIHaloArray edge_sync test is designed with 4 processes only"

function print_array(U,proc)
    if rank == proc
        println()
        println("proc: ", proc)
        for j in size(U, 2):-1:1
            println("j: ", j, "\t", U[:,j])
        end
        istr = join(collect(1:size(U, 1)), "    ")
        println("\t  "*"-"^length(istr))
        println("i :\t  ", istr)
        println()
    end
    MPI.Barrier(comm)
end

function rank_data(ni, nj, rank)
    data = collect(1:ni*nj) * (rank + 1) .+ 10
    reshape(data, ni, nj)
end

function test_edge_sync_2darray_2halo_no_periodic()
    topology = CartesianTopology([4,2], [false, false])
    # topology = CartesianTopology([4,2], [true, true])

    nhalo = 2
    ni = 6
    nj = 5

    data = rank_data(ni, nj, rank)
    A = MPIHaloArray(data, topology, nhalo)
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

    updatehalo!(A)

    if rank == 0
        ihi_neighbor_ilo_edge = rank_data(ni,nj,1)[1:nhalo,:]
        jhi_neighbor_jlo_edge = rank_data(ni,nj,4)[:,1:nhalo]
        @test all(A_ihi_halo .== ihi_neighbor_ilo_edge)
        @test all(A_jhi_halo .== jhi_neighbor_jlo_edge)
        @test all(A_ilo_halo .== -1)
        @test all(A_jlo_halo .== -1)
    end
end

function test_edge_sync_2darray_2halo_all_periodic()
    topology = CartesianTopology([4,2], [true, true])

    nhalo = 2
    ni = 6
    nj = 5

    data = rank_data(ni, nj, rank)
    A = MPIHaloArray(data, topology, nhalo)
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

    updatehalo!(A)
  
    if rank == 0
        ihi_neighbor_ilo_edge = rank_data(ni,nj,1)[1:nhalo,:]
        jhi_neighbor_jlo_edge = rank_data(ni,nj,4)[:,1:nhalo]
        ilo_neighbor_ihi_edge = rank_data(ni,nj,3)[end-nhalo+1:end,:]
        jlo_neighbor_jhi_edge = rank_data(ni,nj,4)[:,end-nhalo+1:end]

        @test all(A_ihi_halo .== ihi_neighbor_ilo_edge)
        @test all(A_jhi_halo .== jhi_neighbor_jlo_edge)
        @test all(A_jlo_halo .== jlo_neighbor_jhi_edge)
        @test all(A_ilo_halo .== ilo_neighbor_ihi_edge)
    end
end

test_edge_sync_2darray_2halo_no_periodic()
test_edge_sync_2darray_2halo_all_periodic()

GC.gc()
MPI.Finalize()
