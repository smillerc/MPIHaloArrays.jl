using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)


@assert nprocs == 4 "This MPIHaloArray edge_sync test is designed with 4 processes only"

function test_edge_sync_2darray_2halo_no_periodic()
    topology = CartesianTopology([2,2], [false, false])

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

    # if rank == 1
    #     println("Before exchange")
    #     @show A.topology.coords
    #     display(A)

    #     println()
    # end
    # println("After exchange")
    # display(A)
    # println()
    # display(A[1,:])

    sync_edges!(A)

    if rank == 0
        @test all(A_ilo_halo .== -1)
        @test all(A_ihi_halo .== 1)
        @test all(A_jlo_halo .== -1)
        @test all(A_jhi_halo .== 2)
    elseif rank == 1
        @test all(A_ilo_halo .== 0)
        @test all(A_ihi_halo .== -1)
        @test all(A_jlo_halo .== -1)
        @test all(A_jhi_halo .== 3)
    elseif rank == 2
        @test all(A_ilo_halo .== -1)
        @test all(A_ihi_halo .== 3)
        @test all(A_jlo_halo .== 0)
        @test all(A_jhi_halo .== -1)
    elseif rank == 3
        @test all(A_ilo_halo .== 2)
        @test all(A_ihi_halo .== -1)
        @test all(A_jlo_halo .== 1)
        @test all(A_jhi_halo .== -1)
    end
end

# function test_edge_sync_2darray_2halo_all_periodic()
#     topology = CartesianTopology([2,2], [true, true])

#     nhalo = 2
#     ni = 8
#     nj = 8


#     A = MPIHaloArray(zeros(Int, ni,nj), topology, nhalo)
#     # Fill each by the current rank
#     filldomain!(A, A.topology.rank)
#     fillhalo!(A, -1)

#     ilo, ihi  = A.local_indices[1].domain
#     jlo, jhi  = A.local_indices[2].domain

#     ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
#     ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
#     jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
#     jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo

#     A_ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi]
#     A_ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi]
#     A_jlo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end]
#     A_jhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end]

#     # if rank == 1
#     #     println("Before exchange")
#     #     @show A.topology.coords
#     #     display(A)

#     #     println()
#     # end
#     # println("After exchange")
#     # display(A)
#     # println()
#     # display(A[1,:])

#     sync_edges!(A)

#     if rank == 0
#         @test all(A_ilo_halo .== -1)
#         @test all(A_ihi_halo .== 1)
#         @test all(A_jlo_halo .== -1)
#         @test all(A_jhi_halo .== 2)
#     elseif rank == 1
#         @test all(A_ilo_halo .== 0)
#         @test all(A_ihi_halo .== -1)
#         @test all(A_jlo_halo .== -1)
#         @test all(A_jhi_halo .== 3)
#     elseif rank == 2
#         @test all(A_ilo_halo .== -1)
#         @test all(A_ihi_halo .== 3)
#         @test all(A_jlo_halo .== 0)
#         @test all(A_jhi_halo .== -1)
#     elseif rank == 3
#         @test all(A_ilo_halo .== 2)
#         @test all(A_ihi_halo .== -1)
#         @test all(A_jlo_halo .== 1)
#         @test all(A_jhi_halo .== -1)
#     end
# end


test_edge_sync_2darray_2halo_no_periodic()


GC.gc()
MPI.Finalize()

# function halo_test()

#     """Helper functions to get the low side halo and domain starting/ending indices"""
#     @inline function lo_indices(nhalo)
#         lo_dom_start = nhalo + 1 # start index of the inner domain
#         lo_dom_end = 2 * nhalo   # end index of the inner domain
#         lo_halo_start = 1        # start index of the halo region
#         lo_halo_end = lo_halo_start + nhalo - 1 # end index of the halo region

#         return (lo_halo_start, lo_halo_end, lo_dom_start, lo_dom_end)
#     end

#     """Helper functions to get the high side halo and domain starting/ending indices"""
#     @inline function hi_indices(field, dim, nhalo)

#         if dim > length(size(field)) error("dim > length(size(field))") end
#         hi_halo_end = size(field, dim)       # end index of the halo region
#         hi_halo_start = hi_halo_end - nhalo + 1  # start index of the halo region

#         hi_dom_end = hi_halo_start - 1   # end index of the inner domain
#         hi_dom_start = hi_dom_end - nhalo + 1 # start index of the inner domain

#         return (hi_dom_start, hi_dom_end, hi_halo_start, hi_halo_end)
#     end

#     function print_array(U)
#         for proc in 0:nprocs
#             if me == proc
#                 println()
#                 println("proc: ", proc)
#                 for j in size(U, 2):-1:1
#                     println("j: ", j, "\t", U[:,j])
#                 end
#                 istr = join(collect(1:size(U, 1)), "   ")
#                 println("\t -"*"-"^length(istr))
#                 println("i :\t  ", istr)
#                 println()
#             end
#             MPI.Barrier(comm)
#         end
#     end

#     ni = 8
#     nj = 10
#     nhalo = 2
#     # U = zeros(ni, nj)
#     # U .= me

#     U = reshape(collect(1:ni*nj) .+ 10, ni, nj) .- 1

#     println("Before")
#     print_array(U)


#     function update_halo_edges!(U)
#         ilo_halo_start, ilo_halo_end, ilo_dom_start, ilo_dom_end = lo_indices(nhalo)
#         jlo_halo_start, jlo_halo_end, jlo_dom_start, jlo_dom_end = lo_indices(nhalo)
#         ihi_dom_start, ihi_dom_end, ihi_halo_start, ihi_halo_end = hi_indices(U, 1, nhalo)
#         jhi_dom_start, jhi_dom_end, jhi_halo_start, jhi_halo_end = hi_indices(U, 2, nhalo)


#         ilo_edge = @view U[ilo_dom_start:ilo_dom_start+nhalo-1, jlo_dom_start:jhi_dom_end]
#         ihi_edge = @view U[ihi_dom_end-nhalo+1:ihi_dom_end, jlo_dom_start:jhi_dom_end]
#         jlo_edge = @view U[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_start+nhalo-1]
#         jhi_edge = @view U[ilo_dom_start:ihi_dom_end, jhi_dom_end-nhalo+1:jhi_dom_end]

#         ilo_halo_edge = @view U[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
#         jlo_halo_edge = @view U[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]

#         ihi_edge_pos = LinearIndices(U)[ihi_dom_end-nhalo+1, jlo_dom_start]
#         ilo_edge_pos = LinearIndices(U)[ilo_dom_start, jlo_dom_start]
#         jlo_edge_pos = LinearIndices(U)[ilo_dom_start, jlo_dom_start]
#         jhi_edge_pos = LinearIndices(U)[ilo_dom_start, jhi_dom_end-nhalo+1]

#         ilo_halo_pos = LinearIndices(U)[ilo_halo_start, jlo_dom_start]
#         ihi_halo_pos = LinearIndices(U)[ihi_halo_start, jlo_dom_start]
#         jlo_halo_pos = LinearIndices(U)[ilo_dom_start, jlo_halo_start]
#         jhi_halo_pos = LinearIndices(U)[ilo_dom_start, jhi_halo_start]

#         # ihi_buf = MPI.Buffer(ihi_edge)
#         # jhi_buf = MPI.Buffer(jhi_edge)
#         ilo_buf = MPI.Buffer(ilo_edge)
#         jlo_buf = MPI.Buffer(jlo_edge)
#         ilo_halo_buf = MPI.Buffer(ilo_halo_edge)
#         jlo_halo_buf = MPI.Buffer(jlo_halo_edge)

#         target_proc = 0

#         ilo_to_ihi_halo = ihi_halo_pos - ilo_edge_pos
#         jlo_to_jhi_halo = jhi_halo_pos - jlo_edge_pos

#         ihi_halo_to_ilo = ihi_halo_pos - ilo_edge_pos
#         jhi_halo_to_jlo = jhi_halo_pos - jlo_edge_pos

#         # Halo exchange
#         win = MPI.Win_create(U, comm)
#         MPI.Win_fence(0, win)
#         MPI.Get(ilo_halo_buf, 1, ihi_halo_to_ilo, win) # update ilo halo; working
#         MPI.Get(jlo_halo_buf, 1, jhi_halo_to_jlo, win) # update the jlo halo; working
#         MPI.Put(jlo_buf, target_proc, jlo_to_jhi_halo, win) # update jhi halo; working
#         MPI.Put(ilo_buf, target_proc, ilo_to_ihi_halo, win) # update ihi halo; working
#         MPI.Win_fence(0, win)
#     end

#     update_halo_edges!(U)

#     ilo_halo = @view U[1:2,3:8]
#     jlo_halo = @view U[3:6,1:2] 
#     ihi_halo = @view U[7:8,3:8]
#     jhi_halo = @view U[3:6,9:10] 

#     ilo_dom = @view U[3:4,3:8]
#     jlo_dom = @view U[3:6,3:4]
#     ihi_dom = @view U[5:6,3:8]
#     jhi_dom = @view U[3:6,7:8]

#     if me == 0
#         @test jhi_halo == jlo_dom
#     elseif me == 1
#         @test jlo_halo == jhi_dom
#     elseif me == 2
#         @test 0 == 1
#     elseif me == 3   
#         @test 0 == 1 
#     end

#     println("After")
#     print_array(U)

#     # MPI.Barrier(MPI.COMM_WORLD)
#     GC.gc()
#     MPI.Finalize()
# end

# halo_test()

