
"""Sync the edges of the array `U` with it's neighbors"""
function sync_edges!(U)
    ilo_halo_start, ilo_halo_end, ilo_dom_start, ilo_dom_end = lo_indices(nhalo)
    jlo_halo_start, jlo_halo_end, jlo_dom_start, jlo_dom_end = lo_indices(nhalo)
    ihi_dom_start, ihi_dom_end, ihi_halo_start, ihi_halo_end = hi_indices(U, 1, nhalo)
    jhi_dom_start, jhi_dom_end, jhi_halo_start, jhi_halo_end = hi_indices(U, 2, nhalo)


    ilo_edge = @view U[ilo_dom_start:ilo_dom_start+nhalo-1, jlo_dom_start:jhi_dom_end]
    ihi_edge = @view U[ihi_dom_end-nhalo+1:ihi_dom_end, jlo_dom_start:jhi_dom_end]
    jlo_edge = @view U[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_start+nhalo-1]
    jhi_edge = @view U[ilo_dom_start:ihi_dom_end, jhi_dom_end-nhalo+1:jhi_dom_end]

    ilo_halo_edge = @view U[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
    jlo_halo_edge = @view U[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]

    ihi_edge_pos = LinearIndices(U)[ihi_dom_end-nhalo+1, jlo_dom_start]
    ilo_edge_pos = LinearIndices(U)[ilo_dom_start, jlo_dom_start]
    jlo_edge_pos = LinearIndices(U)[ilo_dom_start, jlo_dom_start]
    jhi_edge_pos = LinearIndices(U)[ilo_dom_start, jhi_dom_end-nhalo+1]

    ilo_halo_pos = LinearIndices(U)[ilo_halo_start, jlo_dom_start]
    ihi_halo_pos = LinearIndices(U)[ihi_halo_start, jlo_dom_start]
    jlo_halo_pos = LinearIndices(U)[ilo_dom_start, jlo_halo_start]
    jhi_halo_pos = LinearIndices(U)[ilo_dom_start, jhi_halo_start]

    # ihi_buf = MPI.Buffer(ihi_edge)
    # jhi_buf = MPI.Buffer(jhi_edge)
    ilo_buf = MPI.Buffer(ilo_edge)
    jlo_buf = MPI.Buffer(jlo_edge)
    ilo_halo_buf = MPI.Buffer(ilo_halo_edge)
    jlo_halo_buf = MPI.Buffer(jlo_halo_edge)

    target_proc = 0

    ilo_to_ihi_halo = ihi_halo_pos - ilo_edge_pos
    jlo_to_jhi_halo = jhi_halo_pos - jlo_edge_pos

    ihi_halo_to_ilo = ihi_halo_pos - ilo_edge_pos
    jhi_halo_to_jlo = jhi_halo_pos - jlo_edge_pos

    # Halo exchange
    win = MPI.Win_create(U, comm)
    MPI.Win_fence(0, win)
    MPI.Get(ilo_halo_buf, 1, ihi_halo_to_ilo, win) # update ilo halo; working
    MPI.Get(jlo_halo_buf, 1, jhi_halo_to_jlo, win) # update the jlo halo; working
    MPI.Put(jlo_buf, target_proc, jlo_to_jhi_halo, win) # update jhi halo; working
    MPI.Put(ilo_buf, target_proc, ilo_to_ihi_halo, win) # update ihi halo; working
    MPI.Win_fence(0, win)
end
