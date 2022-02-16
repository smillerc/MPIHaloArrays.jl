
# using ..ParallelTopologies

# """Sync the edges of the array `A` with it's neighbors"""
# function sync_edges!(A::MPIHaloArray{T,2}) where {T}

#     nhalo = A.nhalo
#     ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
#     ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
#     ilo_dom_start , ilo_dom_end  = A.local_indices[1].lo_halo_domain_donor
#     ihi_dom_start , ihi_dom_end  = A.local_indices[1].hi_halo_domain_donor

#     jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
#     jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
#     jlo_dom_start , jlo_dom_end  = A.local_indices[2].lo_halo_domain_donor
#     jhi_dom_start , jhi_dom_end  = A.local_indices[2].hi_halo_domain_donor
    

#     # ilo_edge = @view A[ilo_dom_start:ilo_dom_start+nhalo-1, jlo_dom_start:jhi_dom_end]
#     # ihi_edge = @view A[ihi_dom_end-nhalo+1:ihi_dom_end, jlo_dom_start:jhi_dom_end]
#     jlo_edge = @view A[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_start+nhalo-1]
#     # jhi_edge = @view A[ilo_dom_start:ihi_dom_end, jhi_dom_end-nhalo+1:jhi_dom_end]

#     # ilo_halo_edge = @view A[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
#     jlo_halo_edge = @view A[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]

#     # ihi_edge_pos = LinearIndices(A)[ihi_dom_end-nhalo+1, jlo_dom_start]
#     # ilo_edge_pos = LinearIndices(A)[ilo_dom_start, jlo_dom_start]
#     jlo_edge_pos = LinearIndices(A)[ilo_dom_start, jlo_dom_start]
#     jhi_edge_pos = LinearIndices(A)[ilo_dom_start, jhi_dom_end-nhalo+1]

#     # ilo_halo_pos = LinearIndices(A)[ilo_halo_start, jlo_dom_start]
#     # ihi_halo_pos = LinearIndices(A)[ihi_halo_start, jlo_dom_start]
#     jlo_halo_pos = LinearIndices(A)[ilo_dom_start, jlo_halo_start]
#     jhi_halo_pos = LinearIndices(A)[ilo_dom_start, jhi_halo_start]

#     # ihi_buf = MPI.Buffer(ihi_edge)
#     # jhi_buf = MPI.Buffer(jhi_edge)
#     # ilo_buf = MPI.Buffer(ilo_edge)
#     jlo_buf = MPI.Buffer(jlo_edge)
#     # ilo_halo_buf = MPI.Buffer(ilo_halo_edge)
#     jlo_halo_buf = MPI.Buffer(jlo_halo_edge)

#     ilo_neighbor = ilo_neighbor(A.topology)
#     ihi_neighbor = ihi_neighbor(A.topology)
#     jlo_neighbor = jlo_neighbor(A.topology)
#     jhi_neighbor = jhi_neighbor(A.topology)

#     ilo_to_ihi_halo = ihi_halo_pos - ilo_edge_pos
#     jlo_to_jhi_halo = jhi_halo_pos - jlo_edge_pos

#     ihi_halo_to_ilo = ihi_halo_pos - ilo_edge_pos
#     jhi_halo_to_jlo = jhi_halo_pos - jlo_edge_pos

#     # Halo exchange
#     win = MPI.Win_create(A, A.topology.comm)
#     MPI.Win_fence(0, win)
#     # MPI.Get(ilo_halo_buf, target_proc, ihi_halo_to_ilo, win) # update ilo halo; working
#     # MPI.Get(jlo_halo_buf, target_proc, jhi_halo_to_jlo, win) # update the jlo halo; working
#     # MPI.Put(jlo_buf, target_proc, jlo_to_jhi_halo, win) # update jhi halo; working
#     # MPI.Put(ilo_buf, target_proc, ilo_to_ihi_halo, win) # update ihi halo; working
#     MPI.Win_fence(0, win)
# end
