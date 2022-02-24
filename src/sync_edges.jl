"""Sync the edges of the array `A` with it's neighbors"""
function sync_edges_rma!(A::MPIHaloArray{T,1}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start, ilo_dom_end = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start, ihi_dom_end = A.local_indices[1].hi_halo_domain_donor

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end]
    # ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end]

    # Define the positions w/in the window for Get/Put operations
    LI = LinearIndices(A.data)
    ilo_halo_pos = LI[ilo_halo_start]
    ihi_halo_pos = LI[ihi_halo_start]

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

function sync_edges_rma!(A::MPIHaloArray{T,2}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start, ilo_dom_end = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start, ihi_dom_end = A.local_indices[1].hi_halo_domain_donor

    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    jlo_dom_start, jlo_dom_end = A.local_indices[2].lo_halo_domain_donor
    jhi_dom_start, jhi_dom_end = A.local_indices[2].hi_halo_domain_donor

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end]
    jlo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end]

    ilojlo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end]
    ilojhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end]

    ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
    ihi_halo_edge = @view A.data[ihi_halo_start:ihi_halo_end, jlo_dom_start:jhi_dom_end]
    jlo_halo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]

    ilojlo_halo_corner = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end]
    ilojhi_halo_corner = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end]

    # Define the positions w/in the window for Get/Put operations
    LI = LinearIndices(A.data)
    ilo_edge_pos = LI[ilo_dom_start, jlo_dom_start]
    jlo_edge_pos = LI[ilo_dom_start, jlo_dom_start]
    ilojlo_corner_pos = LI[ilo_dom_start, jlo_dom_start]
    ilojhi_corner_pos = LI[ilo_dom_start, jhi_dom_start]

    ihi_halo_pos = LI[ihi_halo_start, jlo_dom_start]
    jhi_halo_pos = LI[ilo_dom_start, jhi_halo_start]
    ihijhi_halo_pos = LI[ihi_halo_start, jhi_halo_start]
    ihijlo_halo_pos = LI[ihi_halo_start, jlo_halo_start]

    ihijhi_corner_pos = LI[ihi_dom_start, jhi_dom_start]

    # Create the MPI subarray buffers for transfer
    ilo_buf = MPI.Buffer(ilo_edge)
    ihi_buf = MPI.Buffer(ihi_edge)
    jlo_buf = MPI.Buffer(jlo_edge)
    ilojlo_corner_buf = MPI.Buffer(ilojlo_corner)
    ilojhi_corner_buf = MPI.Buffer(ilojhi_corner)

    ilo_halo_buf = MPI.Buffer(ilo_halo_edge)
    ihi_halo_buf = MPI.Buffer(ihi_halo_edge)
    jlo_halo_buf = MPI.Buffer(jlo_halo_edge)
    
    ilojlo_halo_buf = MPI.Buffer(ilojlo_halo_corner)
    ilojhi_halo_buf = MPI.Buffer(ilojhi_halo_corner)

    # Calculate the offset to move the buffers
    ilo_to_ihi_halo = ihi_halo_pos - ilo_edge_pos
    ihi_halo_to_ilo = ihi_halo_pos - ilo_edge_pos

    ihijhi_halo_to_ilo_jlo_corner = ihijhi_halo_pos - ilojlo_corner_pos
    ihijlo_halo_to_ilo_jhi_corner = ihijlo_halo_pos - ilojhi_corner_pos

    ihijhi_to_ilojlo_halo = ihijhi_corner_pos - ilojlo_corner_pos

    jlo_to_jhi_halo = jhi_halo_pos - jlo_edge_pos
    jhi_halo_to_jlo = jhi_halo_pos - jlo_edge_pos
    # Halo exchange
    ilo_neighbor_proc = ilo_neighbor(A.topology)
    ihi_neighbor_proc = ihi_neighbor(A.topology)
    jlo_neighbor_proc = jlo_neighbor(A.topology)


    ilojlo_neighbor_proc = neighbor(A.topology,-1,-1)
    ilojhi_neighbor_proc = neighbor(A.topology,-1,+1)

    win = MPI.Win_create(A.data, A.topology.comm)
    MPI.Win_fence(0, win)

    # MPI.Get(ihi_halo_buf, ihi_neighbor_proc, 2, win) # ilo halo update
    # lo halo update: Get the hi edge from the lo neighbor and put it in the lo halo edge

    # Get the hi halo edge from the lo neighbor and put it in the lo halo edge of the current rank
    MPI.Get(ilo_halo_buf, ilo_neighbor_proc, ihi_halo_to_ilo, win) # ilo halo update
    MPI.Get(jlo_halo_buf, jlo_neighbor_proc, jhi_halo_to_jlo, win) # jlo halo update
    MPI.Get(ilojlo_halo_buf, ilojlo_neighbor_proc, ihijhi_to_ilojlo_halo, win) # ilojlo corner halo update
    
    # Put the lo edge (on current rank) on the halo edge of the hi rank
    # MPI.Put(ilo_buf, ilo_neighbor_proc, ilo_to_ihi_halo, win) # ihi halo update
    # MPI.Put(jlo_buf, jlo_neighbor_proc, jlo_to_jhi_halo, win) # jhi halo update
    
    
    # MPI.Get(ilojhi_halo_buf, ilojhi_neighbor_proc, ihijlo_halo_to_ilo_jhi_corner, win) # ilojhi corner halo update
    # MPI.Put(ilojlo_corner_buf, ilojlo_neighbor_proc, jlo_to_jhi_halo, win) # jhi corner halo update
    # MPI.Put(ilojhi_corner_buf, ilojhi_neighbor_proc, ilo_to_ihi_halo, win) # ihi corner halo update

    MPI.Win_fence(0, win)

    return nothing
end

function sync_edges_rma!(A::MPIHaloArray{T,3}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start, ilo_dom_end = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start, ihi_dom_end = A.local_indices[1].hi_halo_domain_donor

    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    jlo_dom_start, jlo_dom_end = A.local_indices[2].lo_halo_domain_donor
    jhi_dom_start, jhi_dom_end = A.local_indices[2].hi_halo_domain_donor

    klo_halo_start, klo_halo_end = A.local_indices[3].lo_halo
    khi_halo_start, khi_halo_end = A.local_indices[3].hi_halo
    klo_dom_start, klo_dom_end = A.local_indices[3].lo_halo_domain_donor
    khi_dom_start, khi_dom_end = A.local_indices[3].hi_halo_domain_donor

    # Edge donor views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end, klo_dom_start:khi_dom_end]
    jlo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:khi_dom_end]
    klo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end]

    # Create the halo region views
    ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end, klo_dom_start:khi_dom_end]
    jlo_halo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end, klo_dom_start:khi_dom_end]
    klo_halo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end, klo_halo_start:klo_halo_end]

    # Define the positions w/in the window for Get/Put operations
    LI = LinearIndices(A.data)
    ilo_edge_pos = LI[ilo_dom_start, jlo_dom_start, klo_dom_start]
    jlo_edge_pos = LI[ilo_dom_start, jlo_dom_start, klo_dom_start]
    klo_edge_pos = LI[ilo_dom_start, jlo_dom_start, klo_dom_start]

    ihi_halo_pos = LI[ihi_halo_start, jlo_dom_start, klo_dom_start]
    jhi_halo_pos = LI[ilo_dom_start, jhi_halo_start, klo_dom_start]
    khi_halo_pos = LI[ilo_dom_start, jhi_dom_start, klo_halo_start]

    if A.topology.rank == 0
        @show ihi_halo_pos, ilo_edge_pos
        @show jhi_halo_pos, jlo_edge_pos
        @show khi_halo_pos, klo_edge_pos
    end
    # Create the MPI subarray buffers for transfer
    ilo_buf = MPI.Buffer(ilo_edge)
    jlo_buf = MPI.Buffer(jlo_edge)
    klo_buf = MPI.Buffer(klo_edge)
    ilo_halo_buf = MPI.Buffer(ilo_halo_edge)
    jlo_halo_buf = MPI.Buffer(jlo_halo_edge)
    klo_halo_buf = MPI.Buffer(klo_halo_edge)

    # Calculate the offset to move the buffers
    ilo_ihi_disp = ihi_halo_pos - ilo_edge_pos
    jlo_jhi_disp = jhi_halo_pos - jlo_edge_pos
    klo_to_khi_halo = khi_halo_pos - klo_edge_pos
    khi_halo_to_klo = khi_halo_pos - klo_edge_pos

    # Halo exchange
    ilo_neighbor_proc = ilo_neighbor(A.topology)
    jlo_neighbor_proc = jlo_neighbor(A.topology)
    klo_neighbor_proc = klo_neighbor(A.topology)

    win = MPI.Win_create(A.data, A.topology.comm)
    MPI.Win_fence(0, win)

    # lo halo update: Get the hi edge from the lo neighbor and put it in the lo halo edge
    MPI.Get(ilo_halo_buf, ilo_neighbor_proc, ilo_ihi_disp, win) # ilo
    MPI.Get(jlo_halo_buf, jlo_neighbor_proc, jlo_jhi_disp, win) # jlo
    # MPI.Get(klo_halo_buf, klo_neighbor_proc, khi_halo_to_klo, win) # klo

    # hi halo updates: Put the lo edge into the lo neighbor's hi halo edge
    MPI.Put(ilo_buf, ilo_neighbor_proc, ilo_ihi_disp, win) # ihi
    MPI.Put(jlo_buf, jlo_neighbor_proc, jlo_jhi_disp, win) # jhi
    # MPI.Put(klo_buf, klo_neighbor_proc, klo_to_khi_halo, win) # khi

    MPI.Win_fence(0, win)

    return nothing
end

function sync_edges!(A::MPIHaloArray{T,1}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start, ilo_dom_end = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start, ihi_dom_end = A.local_indices[1].hi_halo_domain_donor

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end]

    ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end]
    ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end]

    ilo_edge_buf = MPI.Buffer(ilo_edge)
    ihi_edge_buf = MPI.Buffer(ihi_edge)

    ilo_halo_buf = MPI.Buffer(ilo_halo)
    ihi_halo_buf = MPI.Buffer(ihi_halo)

    ilo_neighbor_proc = ilo_neighbor(A.topology)
    ihi_neighbor_proc = ihi_neighbor(A.topology)

    comm = A.topology.comm

    ihi_rreq = MPI.Irecv!(ihi_halo_buf, ihi_neighbor_proc, 1001, comm)
    ilo_rreq = MPI.Irecv!(ilo_halo_buf, ilo_neighbor_proc, 1002, comm)
    ihi_sreq = MPI.Isend(ilo_edge_buf, ilo_neighbor_proc, 1001, comm)
    ilo_sreq = MPI.Isend(ihi_edge_buf, ihi_neighbor_proc, 1002, comm)

    stats = MPI.Waitall!([ihi_rreq, ihi_sreq, ilo_rreq, ilo_sreq])

    return nothing
end

function sync_edges!(A::MPIHaloArray{T,2}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start, ilo_dom_end = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start, ihi_dom_end = A.local_indices[1].hi_halo_domain_donor

    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    jlo_dom_start, jlo_dom_end = A.local_indices[2].lo_halo_domain_donor
    jhi_dom_start, jhi_dom_end = A.local_indices[2].hi_halo_domain_donor

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end]
    jlo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end]
    jhi_edge = @view A.data[ilo_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end]

    ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
    ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_dom_start:jhi_dom_end]
    jlo_halo = @view A.data[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]
    jhi_halo = @view A.data[ilo_dom_start:ihi_dom_end, jhi_halo_start:jhi_halo_end]
   
    ilo_edge_buf = MPI.Buffer(ilo_edge)
    ihi_edge_buf = MPI.Buffer(ihi_edge)
    jlo_edge_buf = MPI.Buffer(jlo_edge)
    jhi_edge_buf = MPI.Buffer(jhi_edge)
    
    ilo_halo_buf = MPI.Buffer(ilo_halo)
    ihi_halo_buf = MPI.Buffer(ihi_halo)
    jlo_halo_buf = MPI.Buffer(jlo_halo)
    jhi_halo_buf = MPI.Buffer(jhi_halo)

    ilo_neighbor_proc = ilo_neighbor(A.topology)
    ihi_neighbor_proc = ihi_neighbor(A.topology)
    jlo_neighbor_proc = jlo_neighbor(A.topology)
    jhi_neighbor_proc = jhi_neighbor(A.topology)

    comm = A.topology.comm

    ihi_rreq = MPI.Irecv!(ihi_halo_buf, ihi_neighbor_proc, 1001, comm)
    ilo_rreq = MPI.Irecv!(ilo_halo_buf, ilo_neighbor_proc, 1002, comm)
    jhi_rreq = MPI.Irecv!(jhi_halo_buf, jhi_neighbor_proc, 1003, comm)
    jlo_rreq = MPI.Irecv!(jlo_halo_buf, jlo_neighbor_proc, 1004, comm)
    
    ihi_sreq = MPI.Isend(ilo_edge_buf, ilo_neighbor_proc, 1001, comm)
    ilo_sreq = MPI.Isend(ihi_edge_buf, ihi_neighbor_proc, 1002, comm)
    jhi_sreq = MPI.Isend(jlo_edge_buf, jlo_neighbor_proc, 1003, comm)
    jlo_sreq = MPI.Isend(jhi_edge_buf, jhi_neighbor_proc, 1004, comm)

    if A.do_corners
        ilojlo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end]
        ilojhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end]
        ihijlo_corner = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end]
        ihijhi_corner = @view A.data[ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end]
        ilojlo_halo   = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end]
        ihijlo_halo   = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end]
        ilojhi_halo   = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end]
        ihijhi_halo   = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end]

        ilojlo_corner_buf = MPI.Buffer(ilojlo_corner)
        ilojhi_corner_buf = MPI.Buffer(ilojhi_corner)
        ihijlo_corner_buf = MPI.Buffer(ihijlo_corner)
        ihijhi_corner_buf = MPI.Buffer(ihijhi_corner)
        ilojlo_halo_buf = MPI.Buffer(ilojlo_halo)  
        ihijlo_halo_buf = MPI.Buffer(ihijlo_halo)  
        ilojhi_halo_buf = MPI.Buffer(ilojhi_halo)  
        ihijhi_halo_buf = MPI.Buffer(ihijhi_halo)  

        ilojlo_neighbor_proc = neighbor(A.topology, -1, -1)
        ihijlo_neighbor_proc = neighbor(A.topology, +1, -1)
        ilojhi_neighbor_proc = neighbor(A.topology, -1, +1)
        ihijhi_neighbor_proc = neighbor(A.topology, +1, +1)

        ilojlo_rreq = MPI.Irecv!(ilojlo_halo_buf  , ilojlo_neighbor_proc, 1005, comm)
        ihijhi_rreq = MPI.Irecv!(ihijhi_halo_buf  , ihijhi_neighbor_proc, 1006, comm)
        ihijlo_rreq = MPI.Irecv!(ihijlo_halo_buf  , ihijlo_neighbor_proc, 1007, comm)
        ilojhi_rreq = MPI.Irecv!(ilojhi_halo_buf  , ilojhi_neighbor_proc, 1008, comm)
        
        ihijhi_sreq =  MPI.Isend(ilojlo_corner_buf, ilojlo_neighbor_proc, 1006, comm)
        ihijlo_sreq =  MPI.Isend(ilojhi_corner_buf, ilojhi_neighbor_proc, 1007, comm)
        ilojlo_sreq =  MPI.Isend(ihijhi_corner_buf, ihijhi_neighbor_proc, 1005, comm)
        ilojhi_sreq =  MPI.Isend(ihijlo_corner_buf, ihijlo_neighbor_proc, 1008, comm)
    end

    if A.do_corners
        stats = MPI.Waitall!([ihi_rreq, ihi_sreq, ilo_rreq, ilo_sreq,
                              jhi_rreq, jhi_sreq, jlo_rreq, jlo_sreq,
                              ilojlo_rreq, ihijlo_rreq, ilojhi_rreq, ihijhi_rreq,
                              ilojlo_sreq, ihijlo_sreq, ilojhi_sreq, ihijhi_sreq])
    else
        stats = MPI.Waitall!([ihi_rreq, ihi_sreq, ilo_rreq, ilo_sreq,
                              jhi_rreq, jhi_sreq, jlo_rreq, jlo_sreq])
    end

    return nothing
end

function sync_edges!(A::MPIHaloArray{T,3}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start, ilo_dom_end = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start, ihi_dom_end = A.local_indices[1].hi_halo_domain_donor

    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    jlo_dom_start, jlo_dom_end = A.local_indices[2].lo_halo_domain_donor
    jhi_dom_start, jhi_dom_end = A.local_indices[2].hi_halo_domain_donor

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end]
    jlo_edge = @view A.data[ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end]
    jhi_edge = @view A.data[ilo_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end]

    ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_dom_start:jhi_dom_end]
    ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_dom_start:jhi_dom_end]
    jlo_halo = @view A.data[ilo_dom_start:ihi_dom_end, jlo_halo_start:jlo_halo_end]
    jhi_halo = @view A.data[ilo_dom_start:ihi_dom_end, jhi_halo_start:jhi_halo_end]
   
    ilo_edge_buf = MPI.Buffer(ilo_edge)
    ihi_edge_buf = MPI.Buffer(ihi_edge)
    jlo_edge_buf = MPI.Buffer(jlo_edge)
    jhi_edge_buf = MPI.Buffer(jhi_edge)
    
    ilo_halo_buf = MPI.Buffer(ilo_halo)
    ihi_halo_buf = MPI.Buffer(ihi_halo)
    jlo_halo_buf = MPI.Buffer(jlo_halo)
    jhi_halo_buf = MPI.Buffer(jhi_halo)

    ilo_neighbor_proc = ilo_neighbor(A.topology)
    ihi_neighbor_proc = ihi_neighbor(A.topology)
    jlo_neighbor_proc = jlo_neighbor(A.topology)
    jhi_neighbor_proc = jhi_neighbor(A.topology)

    comm = A.topology.comm

    ihi_rreq = MPI.Irecv!(ihi_halo_buf, ihi_neighbor_proc, 1001, comm)
    ilo_rreq = MPI.Irecv!(ilo_halo_buf, ilo_neighbor_proc, 1002, comm)
    jhi_rreq = MPI.Irecv!(jhi_halo_buf, jhi_neighbor_proc, 1003, comm)
    jlo_rreq = MPI.Irecv!(jlo_halo_buf, jlo_neighbor_proc, 1004, comm)
    
    ihi_sreq = MPI.Isend(ilo_edge_buf, ilo_neighbor_proc, 1001, comm)
    ilo_sreq = MPI.Isend(ihi_edge_buf, ihi_neighbor_proc, 1002, comm)
    jhi_sreq = MPI.Isend(jlo_edge_buf, jlo_neighbor_proc, 1003, comm)
    jlo_sreq = MPI.Isend(jhi_edge_buf, jhi_neighbor_proc, 1004, comm)

    if A.do_corners
        ilojlo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end]
        ilojhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end]
        ihijlo_corner = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end]
        ihijhi_corner = @view A.data[ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end]
        ilojlo_halo   = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end]
        ihijlo_halo   = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end]
        ilojhi_halo   = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end]
        ihijhi_halo   = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end]

        ilojlo_corner_buf = MPI.Buffer(ilojlo_corner)
        ilojhi_corner_buf = MPI.Buffer(ilojhi_corner)
        ihijlo_corner_buf = MPI.Buffer(ihijlo_corner)
        ihijhi_corner_buf = MPI.Buffer(ihijhi_corner)
        ilojlo_halo_buf = MPI.Buffer(ilojlo_halo)  
        ihijlo_halo_buf = MPI.Buffer(ihijlo_halo)  
        ilojhi_halo_buf = MPI.Buffer(ilojhi_halo)  
        ihijhi_halo_buf = MPI.Buffer(ihijhi_halo)  

        ilojlo_neighbor_proc = neighbor(A.topology, -1, -1)
        ihijlo_neighbor_proc = neighbor(A.topology, +1, -1)
        ilojhi_neighbor_proc = neighbor(A.topology, -1, +1)
        ihijhi_neighbor_proc = neighbor(A.topology, +1, +1)

        ilojlo_rreq = MPI.Irecv!(ilojlo_halo_buf  , ilojlo_neighbor_proc, 1005, comm)
        ihijhi_rreq = MPI.Irecv!(ihijhi_halo_buf  , ihijhi_neighbor_proc, 1006, comm)
        ihijlo_rreq = MPI.Irecv!(ihijlo_halo_buf  , ihijlo_neighbor_proc, 1007, comm)
        ilojhi_rreq = MPI.Irecv!(ilojhi_halo_buf  , ilojhi_neighbor_proc, 1008, comm)
        
        ihijhi_sreq =  MPI.Isend(ilojlo_corner_buf, ilojlo_neighbor_proc, 1006, comm)
        ihijlo_sreq =  MPI.Isend(ilojhi_corner_buf, ilojhi_neighbor_proc, 1007, comm)
        ilojlo_sreq =  MPI.Isend(ihijhi_corner_buf, ihijhi_neighbor_proc, 1005, comm)
        ilojhi_sreq =  MPI.Isend(ihijlo_corner_buf, ihijlo_neighbor_proc, 1008, comm)
    end

    if A.do_corners
        stats = MPI.Waitall!([ihi_rreq, ihi_sreq, ilo_rreq, ilo_sreq,
                              jhi_rreq, jhi_sreq, jlo_rreq, jlo_sreq,
                              ilojlo_rreq, ihijlo_rreq, ilojhi_rreq, ihijhi_rreq,
                              ilojlo_sreq, ihijlo_sreq, ilojhi_sreq, ihijhi_sreq])
    else
        stats = MPI.Waitall!([ihi_rreq, ihi_sreq, ilo_rreq, ilo_sreq,
                              jhi_rreq, jhi_sreq, jlo_rreq, jlo_sreq])
    end

    return nothing
end