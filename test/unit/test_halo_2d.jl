using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This MPIHaloArray edge_sync test is designed with 4 processes only"

function print_array(U, proc)
    if rank == proc
        println()
        println("proc: ", proc)
        for j in size(U, 2):-1:1
            println("j: ", j, "\t", U[:, j])
        end
        istr = join(collect(1:size(U, 1)), "    ")
        println("\t  " * "-"^length(istr))
        println("i :\t  ", istr)
        println()
    end
    MPI.Barrier(comm)

    return nothing
end

function rank_data(ni, nj, rank)
    data = collect(1:ni*nj) * (rank + 1) .+ 10
    return reshape(data, ni, nj)
end

"""Make a 3D array designed to be used with 2 halo dimensions only"""
function rank_data3d_2d_halo_dims(n, ni, nj, rank)

    data = zeros(Int, n, ni, nj)

    for j in 1:nj
        for i in 1:ni
            v = i * j * (rank + 1) + 10 # same logic as in the `rank_data` function, now with an extra dimension
            for q in 1:n
                data[q, i, j] = v - (2q - 1)
            end
        end
    end

    # data[1,:,:] will have the same results as `rank_data`
    return data
end

function test_edge_sync_2darray_2halo_no_periodic()
    topology = CartesianTopology(comm, (4, 2), (false, false))

    nhalo = 2
    ni = 6
    nj = 5

    data = rank_data(ni, nj, rank)
    A = MPIHaloArray(data, topology, nhalo)
    fillhalo!(A, -1)

    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain

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
        ihi_neighbor_ilo_edge = rank_data(ni, nj, 1)[1:nhalo, :]
        jhi_neighbor_jlo_edge = rank_data(ni, nj, 4)[:, 1:nhalo]
        @test all(A_ihi_halo .== ihi_neighbor_ilo_edge)
        @test all(A_jhi_halo .== jhi_neighbor_jlo_edge)
        @test all(A_ilo_halo .== -1)
        @test all(A_jlo_halo .== -1)
    end
end

function test_edge_sync_3darray_2halodims_2halocells_no_periodic()

    # Rank layout
    #   +---- (i)
    #   |
    #   | (j)
    #
    #  jlo   |  0  |  1  |  2  |  3  |
    #  jhi   |  4  |  5  |  6  |  7  |
    #          ilo               ihi

    topology = CartesianTopology(comm, (4, 2), (false, false))

    nhalo = 2
    ni = 6
    nj = 5
    N = 4
    halo_dims = (2, 3) # only do halo exchange on the 2nd and 3rd dimensions

    data = rank_data3d_2d_halo_dims(N, ni, nj, rank)
    A = MPIHaloArray(data, topology, nhalo, halo_dims)
    fillhalo!(A, -1)

    # ni is 6, so i goes 1:6, but with 2 halo cells, it goes 1 2 | 3 4 5 6 7 8 | 9 10,
    # so ilo should be 3 and ihi should be 8

    # ignore the first dimension, because no halo exchange is done on it
    _, _, ilo, ihi, jlo, jhi = local_domain_indices(A)
    @test (ilo, ihi, jlo, jhi) == (3, 8, 3, 7)

    ilo_halo_start, ilo_halo_end = A.local_indices[2].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[2].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[3].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[3].hi_halo



    A_ilo_halo = @view A.data[:, ilo_halo_start:ilo_halo_end, jlo:jhi]
    A_ihi_halo = @view A.data[:, ihi_halo_start:ihi_halo_end, jlo:jhi]
    A_jlo_halo = @view A.data[:, ilo:ihi, jlo_halo_start:jlo_halo_end]
    A_jhi_halo = @view A.data[:, ilo:ihi, jhi_halo_start:jhi_halo_end]

    updatehalo!(A)

    if rank == 0
        # # Rank 0's neighbor on the ihi side is 1
        ihi_neighbor_rank = 1
        ihi_neighbor_data = rank_data3d_2d_halo_dims(N, ni, nj, ihi_neighbor_rank)
        ihi_neighbor_ilo_edge = ihi_neighbor_data[:, 1:nhalo, :]
        @test all(A_ihi_halo .== ihi_neighbor_ilo_edge)

        # Rank 0's neighbor on the jhi side is 4
        ihi_neighbor = 4
        jhi_neighbor_jlo_data = rank_data3d_2d_halo_dims(4, ni, nj, ihi_neighbor)
        jhi_neighbor_jlo_edge = jhi_neighbor_jlo_data[:, :, 1:nhalo]
        @test all(A_jhi_halo .== jhi_neighbor_jlo_edge)
        @test all(A_ilo_halo .== -1)
        @test all(A_jlo_halo .== -1)
    end
end

function test_edge_sync_3darray_2halodims_2halocells_all_periodic()

    # Rank layout
    #   +---- (i)
    #   |
    #   | (j)
    #
    #  jlo   |  0  |  1  |  2  |  3  |
    #  jhi   |  4  |  5  |  6  |  7  |
    #          ilo               ihi

    topology = CartesianTopology(comm, (4, 2), (true, true))

    nhalo = 2
    ni = 6
    nj = 5
    N = 4
    halo_dims = (2, 3) # only do halo exchange on the 2nd and 3rd dimensions

    data = rank_data3d_2d_halo_dims(N, ni, nj, rank)
    A = MPIHaloArray(data, topology, nhalo, halo_dims)
    fillhalo!(A, -1)

    # ni is 6, so i goes 1:6, but with 2 halo cells, it goes 1 2 | 3 4 5 6 7 8 | 9 10,
    # so ilo should be 3 and ihi should be 8

    # ignore the first dimension, because no halo exchange is done on it
    _, _, ilo, ihi, jlo, jhi = local_domain_indices(A)
    @test (ilo, ihi, jlo, jhi) == (3, 8, 3, 7)

    ilo_halo_start, ilo_halo_end = A.local_indices[2].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[2].hi_halo
    jlo_halo_start, jlo_halo_end = A.local_indices[3].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[3].hi_halo



    A_ilo_halo = @view A.data[:, ilo_halo_start:ilo_halo_end, jlo:jhi]
    A_ihi_halo = @view A.data[:, ihi_halo_start:ihi_halo_end, jlo:jhi]
    A_jlo_halo = @view A.data[:, ilo:ihi, jlo_halo_start:jlo_halo_end]
    A_jhi_halo = @view A.data[:, ilo:ihi, jhi_halo_start:jhi_halo_end]

    updatehalo!(A)

    if rank == 0
        ihi_neighbor_ilo_edge = rank_data3d_2d_halo_dims(N, ni, nj, 1)[:, 1:nhalo, :]
        jhi_neighbor_jlo_edge = rank_data3d_2d_halo_dims(N, ni, nj, 4)[:, :, 1:nhalo]
        ilo_neighbor_ihi_edge = rank_data3d_2d_halo_dims(N, ni, nj, 3)[:, end-nhalo+1:end, :]
        jlo_neighbor_jhi_edge = rank_data3d_2d_halo_dims(N, ni, nj, 4)[:, :, end-nhalo+1:end]

        @test all(A_ihi_halo .== ihi_neighbor_ilo_edge)
        @test all(A_jhi_halo .== jhi_neighbor_jlo_edge)
        @test all(A_jlo_halo .== jlo_neighbor_jhi_edge)
        @test all(A_ilo_halo .== ilo_neighbor_ihi_edge)
    end
end

function test_edge_sync_2darray_2halo_all_periodic()
    topology = CartesianTopology(comm, (4, 2), (true, true))

    nhalo = 2
    ni = 6
    nj = 5

    data = rank_data(ni, nj, rank)
    A = MPIHaloArray(data, topology, nhalo)
    fillhalo!(A, -1)

    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain

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
        ihi_neighbor_ilo_edge = rank_data(ni, nj, 1)[1:nhalo, :]
        jhi_neighbor_jlo_edge = rank_data(ni, nj, 4)[:, 1:nhalo]
        ilo_neighbor_ihi_edge = rank_data(ni, nj, 3)[end-nhalo+1:end, :]
        jlo_neighbor_jhi_edge = rank_data(ni, nj, 4)[:, end-nhalo+1:end]

        @test all(A_ihi_halo .== ihi_neighbor_ilo_edge)
        @test all(A_jhi_halo .== jhi_neighbor_jlo_edge)
        @test all(A_jlo_halo .== jlo_neighbor_jhi_edge)
        @test all(A_ilo_halo .== ilo_neighbor_ihi_edge)
    end
end

test_edge_sync_2darray_2halo_no_periodic()
test_edge_sync_2darray_2halo_all_periodic()
test_edge_sync_3darray_2halodims_2halocells_no_periodic()
test_edge_sync_3darray_2halodims_2halocells_all_periodic()

GC.gc()
MPI.Finalize()
