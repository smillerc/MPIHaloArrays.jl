"""Sync the edges of the array `A` with it's neighbors"""
function sync_edges!(A::MPIHaloArray{T,1}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start , ilo_dom_end  = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start , ihi_dom_end  = A.local_indices[1].hi_halo_domain_donor    

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end]
    # ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end]

    # Define the positions w/in the window for Get/Put operations
    LI = LinearIndices(A.data)
    ilo_halo_pos =LI[ilo_halo_start]
    ihi_halo_pos =LI[ihi_halo_start]

    # Create the MPI subarray buffers for transfer
    ilo_buf = MPI.Buffer(ilo_edge)
    ihi_buf = MPI.Buffer(ihi_edge)

    # Halo exchange
    ilo_neighbor_proc = ilo_neighbor(A.topology)
    ihi_neighbor_proc = ihi_neighbor(A.topology)

    ilo_to_ihi_halo_disp = ihi_halo_pos - ilo_halo_pos
    ihi_to_ilo_halo_disp = 0

    # if A.topology.rank == 1
    #     @show ilo_edge_pos, ihi_edge_pos, ilo_halo_pos, ihi_halo_pos
    #     @show ihi_halo_pos - ilo_dom_start + 2
    # end

    win = MPI.Win_create(A.data, A.topology.comm)
    MPI.Win_fence(0, win)
    MPI.Put(ilo_buf, ilo_neighbor_proc, ilo_to_ihi_halo_disp, win) # ihi halo update
    MPI.Put(ihi_buf, ihi_neighbor_proc, ihi_to_ilo_halo_disp, win) # ilo halo update
    MPI.Win_fence(0, win)

    return nothing
end

function sync_edges!(A::MPIHaloArray{T,2}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start , ilo_dom_end  = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start , ihi_dom_end  = A.local_indices[1].hi_halo_domain_donor

    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    jlo_dom_start , jlo_dom_end  = A.local_indices[2].lo_halo_domain_donor
    jhi_dom_start , jhi_dom_end  = A.local_indices[2].hi_halo_domain_donor
    

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end]
    jlo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end]
    ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
    jlo_halo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]

    # Define the positions w/in the window for Get/Put operations
    LI = LinearIndices(A.data)
    ilo_edge_pos =LI[ilo_dom_start, jlo_dom_start]
    jlo_edge_pos =LI[ilo_dom_start, jlo_dom_start]
    ihi_halo_pos =LI[ihi_halo_start, jlo_dom_start]
    jhi_halo_pos =LI[ilo_dom_start, jhi_halo_start]

    # Create the MPI subarray buffers for transfer
    ilo_buf = MPI.Buffer(ilo_edge)
    jlo_buf = MPI.Buffer(jlo_edge)
    ilo_halo_buf = MPI.Buffer(ilo_halo_edge)
    jlo_halo_buf = MPI.Buffer(jlo_halo_edge)

    # Calculate the offset to move the buffers
    ilo_to_ihi_halo = ihi_halo_pos - ilo_edge_pos
    ihi_halo_to_ilo = ihi_halo_pos - ilo_edge_pos
    jlo_to_jhi_halo = jhi_halo_pos - jlo_edge_pos
    jhi_halo_to_jlo = jhi_halo_pos - jlo_edge_pos

    # Halo exchange
    ilo_neighbor_proc = ilo_neighbor(A.topology)
    jlo_neighbor_proc = jlo_neighbor(A.topology)

    win = MPI.Win_create(A.data, A.topology.comm)
    MPI.Win_fence(0, win)

    MPI.Get(ilo_halo_buf, ilo_neighbor_proc, ihi_halo_to_ilo, win) # ilo halo update
    MPI.Get(jlo_halo_buf, jlo_neighbor_proc, jhi_halo_to_jlo, win) # jlo halo update
    MPI.Put(jlo_buf,      jlo_neighbor_proc, jlo_to_jhi_halo, win) # jhi halo update
    MPI.Put(ilo_buf,      ilo_neighbor_proc, ilo_to_ihi_halo, win) # ihi halo update
    
    MPI.Win_fence(0, win)

    return nothing
end
