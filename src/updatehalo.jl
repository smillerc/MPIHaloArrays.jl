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
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo:jhi]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end, jlo:jhi]
    jlo_edge = @view A.data[ilo:ihi, jlo_dom_start:jlo_dom_end]

    ilojlo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end]
    ilojhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end]

    ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi]
    ihi_halo_edge = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi]
    jlo_halo_edge = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end]

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
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo:jhi, klo:khi]
    jlo_edge = @view A.data[ilo:ihi, jlo_dom_start:jlo_dom_end, klo:khi]
    klo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo:jhi, klo_dom_start:klo_dom_end]

    # Create the halo region views
    ilo_halo_edge = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, klo:khi]
    jlo_halo_edge = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, klo:khi]
    klo_halo_edge = @view A.data[ilo:ihi, jlo:jhi, klo_halo_start:klo_halo_end]

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

function updatehalo!(A::MPIHaloArray{T,1}) where {T}
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

function updatehalo!(A::MPIHaloArray{T,2}) where {T}
    ilo_halo_start, ilo_halo_end = A.local_indices[1].lo_halo
    ihi_halo_start, ihi_halo_end = A.local_indices[1].hi_halo
    ilo_dom_start, ilo_dom_end = A.local_indices[1].lo_halo_domain_donor
    ihi_dom_start, ihi_dom_end = A.local_indices[1].hi_halo_domain_donor

    jlo_halo_start, jlo_halo_end = A.local_indices[2].lo_halo
    jhi_halo_start, jhi_halo_end = A.local_indices[2].hi_halo
    jlo_dom_start, jlo_dom_end = A.local_indices[2].lo_halo_domain_donor
    jhi_dom_start, jhi_dom_end = A.local_indices[2].hi_halo_domain_donor

    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo:jhi]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end, jlo:jhi]
    jlo_edge = @view A.data[ilo:ihi, jlo_dom_start:jlo_dom_end]
    jhi_edge = @view A.data[ilo:ihi, jhi_dom_start:jhi_dom_end]

    ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi]
    ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi]
    jlo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end]
    jhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end]

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

function updatehalo!(A::MPIHaloArray{T,3}) where {T}
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

    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain
    klo, khi = A.local_indices[3].domain

    # Create the halo region views
    ilo_edge = @view A.data[ilo_dom_start:ilo_dom_end, jlo:jhi, klo:khi]
    ihi_edge = @view A.data[ihi_dom_start:ihi_dom_end, jlo:jhi, klo:khi]
    jlo_edge = @view A.data[ilo:ihi, jlo_dom_start:jlo_dom_end, klo:khi]
    jhi_edge = @view A.data[ilo:ihi, jhi_dom_start:jhi_dom_end, klo:khi]
    klo_edge = @view A.data[ilo:ihi, jlo:jhi, klo_dom_start:klo_dom_end]
    khi_edge = @view A.data[ilo:ihi, jlo:jhi, khi_dom_start:khi_dom_end]

    ilo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, klo:khi]
    ihi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, klo:khi]
    jlo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, klo:khi]
    jhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, klo:khi]
    klo_halo = @view A.data[ilo:ihi, jlo:jhi, klo_halo_start:klo_halo_end]
    khi_halo = @view A.data[ilo:ihi, jlo:jhi, khi_halo_start:khi_halo_end]

    ilo_edge_buf = MPI.Buffer(ilo_edge)
    ihi_edge_buf = MPI.Buffer(ihi_edge)
    jlo_edge_buf = MPI.Buffer(jlo_edge)
    jhi_edge_buf = MPI.Buffer(jhi_edge)
    klo_edge_buf = MPI.Buffer(klo_edge)
    khi_edge_buf = MPI.Buffer(khi_edge)

    ilo_halo_buf = MPI.Buffer(ilo_halo)
    ihi_halo_buf = MPI.Buffer(ihi_halo)
    jlo_halo_buf = MPI.Buffer(jlo_halo)
    jhi_halo_buf = MPI.Buffer(jhi_halo)
    klo_halo_buf = MPI.Buffer(klo_halo)
    khi_halo_buf = MPI.Buffer(khi_halo)

    ilo_neighbor_proc = ilo_neighbor(A.topology)
    ihi_neighbor_proc = ihi_neighbor(A.topology)
    jlo_neighbor_proc = jlo_neighbor(A.topology)
    jhi_neighbor_proc = jhi_neighbor(A.topology)
    klo_neighbor_proc = klo_neighbor(A.topology)
    khi_neighbor_proc = khi_neighbor(A.topology)

    comm = A.topology.comm

    ihi_rreq = MPI.Irecv!(ihi_halo_buf, ihi_neighbor_proc, 1001, comm)
    ilo_rreq = MPI.Irecv!(ilo_halo_buf, ilo_neighbor_proc, 1002, comm)
    jhi_rreq = MPI.Irecv!(jhi_halo_buf, jhi_neighbor_proc, 1003, comm)
    jlo_rreq = MPI.Irecv!(jlo_halo_buf, jlo_neighbor_proc, 1004, comm)
    khi_rreq = MPI.Irecv!(khi_halo_buf, khi_neighbor_proc, 1005, comm)
    klo_rreq = MPI.Irecv!(klo_halo_buf, klo_neighbor_proc, 1006, comm)

    ihi_sreq = MPI.Isend(ilo_edge_buf, ilo_neighbor_proc, 1001, comm)
    ilo_sreq = MPI.Isend(ihi_edge_buf, ihi_neighbor_proc, 1002, comm)
    jhi_sreq = MPI.Isend(jlo_edge_buf, jlo_neighbor_proc, 1003, comm)
    jlo_sreq = MPI.Isend(jhi_edge_buf, jhi_neighbor_proc, 1004, comm)
    khi_sreq = MPI.Isend(klo_edge_buf, klo_neighbor_proc, 1005, comm)
    klo_sreq = MPI.Isend(khi_edge_buf, khi_neighbor_proc, 1006, comm)

    if A.do_corners
        jloklo_corner = @view A.data[ilo:ihi, jlo_dom_start:jlo_dom_end, klo_dom_start:klo_dom_end]
        jlokhi_corner = @view A.data[ilo:ihi, jlo_dom_start:jlo_dom_end, khi_dom_start:khi_dom_end]
        jhiklo_corner = @view A.data[ilo:ihi, jhi_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end]
        jhikhi_corner = @view A.data[ilo:ihi, jhi_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end]
        
        iloklo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo:jhi, klo_dom_start:klo_dom_end]
        ilokhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo:jhi, khi_dom_start:khi_dom_end]
        ihiklo_corner = @view A.data[ihi_dom_start:ihi_dom_end, jlo:jhi, klo_dom_start:klo_dom_end]
        ihikhi_corner = @view A.data[ihi_dom_start:ihi_dom_end, jlo:jhi, khi_dom_start:khi_dom_end]
        
        ilojlo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, klo:khi]
        ilojhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, klo:khi]
        ihijlo_corner = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo:khi]
        ihijhi_corner = @view A.data[ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, klo:khi]

        ilojloklo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:klo_dom_end]
        ilojhiklo_corner = @view A.data[ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end]
        ihijloklo_corner = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:klo_dom_end]
        ihijhiklo_corner = @view A.data[ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end]

        ilojlokhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, khi_dom_start:khi_dom_end]
        ilojhikhi_corner = @view A.data[ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end]
        ihijlokhi_corner = @view A.data[ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, khi_dom_start:khi_dom_end]
        ihijhikhi_corner = @view A.data[ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end]

        ilojlo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo:khi]
        ilojhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo:khi]
        ihijlo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo:khi]
        ihijhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo:khi]
        
        jloklo_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
        jlokhi_halo = @view A.data[ilo:ihi, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
        jhiklo_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
        jhikhi_halo = @view A.data[ilo:ihi, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]
        
        iloklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, klo_halo_start:klo_halo_end]
        ilokhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo:jhi, khi_halo_start:khi_halo_end]
        ihiklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, klo_halo_start:klo_halo_end]
        ihikhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo:jhi, khi_halo_start:khi_halo_end]

        ilojloklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
        ilojhiklo_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]
        ihijloklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end]
        ihijhiklo_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end]

        ilojlokhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
        ilojhikhi_halo = @view A.data[ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]
        ihijlokhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end]
        ihijhikhi_halo = @view A.data[ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end]

        ilojlo_corner_buf = MPI.Buffer(ilojlo_corner)
        ilojhi_corner_buf = MPI.Buffer(ilojhi_corner)
        ihijlo_corner_buf = MPI.Buffer(ihijlo_corner)
        ihijhi_corner_buf = MPI.Buffer(ihijhi_corner)

        jloklo_corner_buf = MPI.Buffer(jloklo_corner)
        jlokhi_corner_buf = MPI.Buffer(jlokhi_corner)
        jhiklo_corner_buf = MPI.Buffer(jhiklo_corner)
        jhikhi_corner_buf = MPI.Buffer(jhikhi_corner)
        
        iloklo_corner_buf = MPI.Buffer(iloklo_corner)
        ilokhi_corner_buf = MPI.Buffer(ilokhi_corner)
        ihiklo_corner_buf = MPI.Buffer(ihiklo_corner)
        ihikhi_corner_buf = MPI.Buffer(ihikhi_corner)

        ilojloklo_corner_buf = MPI.Buffer(ilojloklo_corner)
        ilojhiklo_corner_buf = MPI.Buffer(ilojhiklo_corner)
        ihijloklo_corner_buf = MPI.Buffer(ihijloklo_corner)
        ihijhiklo_corner_buf = MPI.Buffer(ihijhiklo_corner)

        ilojlokhi_corner_buf = MPI.Buffer(ilojlokhi_corner)
        ilojhikhi_corner_buf = MPI.Buffer(ilojhikhi_corner)
        ihijlokhi_corner_buf = MPI.Buffer(ihijlokhi_corner)
        ihijhikhi_corner_buf = MPI.Buffer(ihijhikhi_corner)

        ilojlo_halo_buf = MPI.Buffer(ilojlo_halo)
        ilojhi_halo_buf = MPI.Buffer(ilojhi_halo)
        ihijlo_halo_buf = MPI.Buffer(ihijlo_halo)
        ihijhi_halo_buf = MPI.Buffer(ihijhi_halo)

        jloklo_halo_buf = MPI.Buffer(jloklo_halo)
        jlokhi_halo_buf = MPI.Buffer(jlokhi_halo)
        jhiklo_halo_buf = MPI.Buffer(jhiklo_halo)
        jhikhi_halo_buf = MPI.Buffer(jhikhi_halo)
        
        iloklo_halo_buf = MPI.Buffer(iloklo_halo)
        ilokhi_halo_buf = MPI.Buffer(ilokhi_halo)
        ihiklo_halo_buf = MPI.Buffer(ihiklo_halo)
        ihikhi_halo_buf = MPI.Buffer(ihikhi_halo)

        ilojloklo_halo_buf = MPI.Buffer(ilojloklo_halo)
        ilojhiklo_halo_buf = MPI.Buffer(ilojhiklo_halo)
        ihijloklo_halo_buf = MPI.Buffer(ihijloklo_halo)
        ihijhiklo_halo_buf = MPI.Buffer(ihijhiklo_halo)

        ilojlokhi_halo_buf = MPI.Buffer(ilojlokhi_halo)
        ilojhikhi_halo_buf = MPI.Buffer(ilojhikhi_halo)
        ihijlokhi_halo_buf = MPI.Buffer(ihijlokhi_halo)
        ihijhikhi_halo_buf = MPI.Buffer(ihijhikhi_halo)

        ilojlo_neighbor_proc = neighbor(A.topology, -1, -1, 0)
        ihijlo_neighbor_proc = neighbor(A.topology, +1, -1, 0)
        ilojhi_neighbor_proc = neighbor(A.topology, -1, +1, 0)
        ihijhi_neighbor_proc = neighbor(A.topology, +1, +1, 0)
        
        iloklo_neighbor_proc = neighbor(A.topology, -1, 0, -1)
        ihiklo_neighbor_proc = neighbor(A.topology, +1, 0, -1)
        ilokhi_neighbor_proc = neighbor(A.topology, -1, 0, +1)
        ihikhi_neighbor_proc = neighbor(A.topology, +1, 0, +1)
        
        jloklo_neighbor_proc = neighbor(A.topology, 0, -1, -1)
        jhiklo_neighbor_proc = neighbor(A.topology, 0, +1, -1)
        jlokhi_neighbor_proc = neighbor(A.topology, 0, -1, +1)
        jhikhi_neighbor_proc = neighbor(A.topology, 0, +1, +1)

        ilojloklo_neighbor_proc = neighbor(A.topology, -1,-1,-1)
        ilojhiklo_neighbor_proc = neighbor(A.topology, -1,+1,-1)
        ihijloklo_neighbor_proc = neighbor(A.topology, +1,-1,-1)
        ihijhiklo_neighbor_proc = neighbor(A.topology, +1,+1,-1)
        ilojlokhi_neighbor_proc = neighbor(A.topology, -1,-1,+1)
        ilojhikhi_neighbor_proc = neighbor(A.topology, -1,+1,+1)
        ihijlokhi_neighbor_proc = neighbor(A.topology, +1,-1,+1)
        ihijhikhi_neighbor_proc = neighbor(A.topology, +1,+1,+1)

        ilojlo_rreq    = MPI.Irecv!(ilojlo_halo_buf    , ilojlo_neighbor_proc   , 1007, comm)
        ilojlo_sreq    =  MPI.Isend(ihijhi_corner_buf  , ihijhi_neighbor_proc   , 1007, comm)

        ihijhi_rreq    = MPI.Irecv!(ihijhi_halo_buf    , ihijhi_neighbor_proc   , 1008, comm)
        ihijhi_sreq    =  MPI.Isend(ilojlo_corner_buf  , ilojlo_neighbor_proc   , 1008, comm)

        ihijlo_rreq    = MPI.Irecv!(ihijlo_halo_buf    , ihijlo_neighbor_proc   , 1009, comm)
        ihijlo_sreq    =  MPI.Isend(ilojhi_corner_buf  , ilojhi_neighbor_proc   , 1009, comm)

        ilojhi_rreq    = MPI.Irecv!(ilojhi_halo_buf    , ilojhi_neighbor_proc   , 1010, comm)
        ilojhi_sreq    =  MPI.Isend(ihijlo_corner_buf  , ihijlo_neighbor_proc   , 1010, comm)



        iloklo_halo_rreq = MPI.Irecv!(iloklo_halo_buf, iloklo_neighbor_proc, 1011, comm)
        iloklo_halo_sreq =  MPI.Isend(ihikhi_corner_buf, ihikhi_neighbor_proc, 1011, comm)

        ilokhi_halo_rreq = MPI.Irecv!(ilokhi_halo_buf, ilokhi_neighbor_proc, 1012, comm)
        ilokhi_halo_sreq =  MPI.Isend(ihiklo_corner_buf, ihiklo_neighbor_proc, 1012, comm)

        ihiklo_halo_rreq = MPI.Irecv!(ihiklo_halo_buf, ihiklo_neighbor_proc, 1013, comm)
        ihiklo_halo_sreq =  MPI.Isend(ilokhi_corner_buf, ilokhi_neighbor_proc, 1013, comm)

        ihikhi_halo_rreq = MPI.Irecv!(ihikhi_halo_buf, ihikhi_neighbor_proc, 1014, comm)
        ihikhi_halo_sreq =  MPI.Isend(iloklo_corner_buf, iloklo_neighbor_proc, 1014, comm)

        jloklo_halo_rreq = MPI.Irecv!(jloklo_halo_buf, jloklo_neighbor_proc, 1015, comm)
        jloklo_halo_sreq =  MPI.Isend(jhikhi_corner_buf, jhikhi_neighbor_proc, 1015, comm)

        jlokhi_halo_rreq = MPI.Irecv!(jlokhi_halo_buf, jlokhi_neighbor_proc, 1016, comm)
        jlokhi_halo_sreq =  MPI.Isend(jhiklo_corner_buf, jhiklo_neighbor_proc, 1016, comm)

        jhiklo_halo_rreq = MPI.Irecv!(jhiklo_halo_buf, jhiklo_neighbor_proc, 1017, comm)
        jhiklo_halo_sreq =  MPI.Isend(jlokhi_corner_buf, jlokhi_neighbor_proc, 1017, comm)

        jhikhi_halo_rreq = MPI.Irecv!(jhikhi_halo_buf, jhikhi_neighbor_proc, 1018, comm)
        jhikhi_halo_sreq =  MPI.Isend(jloklo_corner_buf, jloklo_neighbor_proc, 1018, comm)




        ilojloklo_rreq = MPI.Irecv!(ilojloklo_halo_buf,   ilojloklo_neighbor_proc, 1019, comm)
        ilojloklo_sreq =  MPI.Isend(ihijhikhi_corner_buf, ihijhikhi_neighbor_proc, 1019, comm)

        ilojhiklo_rreq = MPI.Irecv!(ilojhiklo_halo_buf,   ilojhiklo_neighbor_proc, 1020, comm)
        ilojhiklo_sreq =  MPI.Isend(ihijlokhi_corner_buf, ihijlokhi_neighbor_proc, 1020, comm)

        ihijloklo_rreq = MPI.Irecv!(ihijloklo_halo_buf,   ihijloklo_neighbor_proc, 1021, comm)
        ihijloklo_sreq =  MPI.Isend(ilojhikhi_corner_buf, ilojhikhi_neighbor_proc, 1021, comm)

        ihijhiklo_rreq = MPI.Irecv!(ihijhiklo_halo_buf,   ihijhiklo_neighbor_proc, 1022, comm)
        ihijhiklo_sreq =  MPI.Isend(ilojlokhi_corner_buf, ilojlokhi_neighbor_proc, 1022, comm)

        ilojlokhi_rreq = MPI.Irecv!(ilojlokhi_halo_buf,   ilojlokhi_neighbor_proc, 1023, comm)
        ilojlokhi_sreq =  MPI.Isend(ihijhiklo_corner_buf, ihijhiklo_neighbor_proc, 1023, comm)

        ilojhikhi_rreq = MPI.Irecv!(ilojhikhi_halo_buf,   ilojhikhi_neighbor_proc, 1024, comm)
        ilojhikhi_sreq =  MPI.Isend(ihijloklo_corner_buf, ihijloklo_neighbor_proc, 1024, comm)

        ihijlokhi_rreq = MPI.Irecv!(ihijlokhi_halo_buf,   ihijlokhi_neighbor_proc, 1025, comm)
        ihijlokhi_sreq =  MPI.Isend(ilojhiklo_corner_buf, ilojhiklo_neighbor_proc, 1025, comm)

        ihijhikhi_rreq = MPI.Irecv!(ihijhikhi_halo_buf,   ihijhikhi_neighbor_proc, 1026, comm)
        ihijhikhi_sreq =  MPI.Isend(ilojloklo_corner_buf, ilojloklo_neighbor_proc, 1026, comm)
    end

    if A.do_corners
        stats = MPI.Waitall!([ihi_rreq, ihi_sreq, ilo_rreq, ilo_sreq,
                             jhi_rreq, jhi_sreq, jlo_rreq, jlo_sreq,
                             khi_rreq, khi_sreq, klo_rreq, klo_sreq,
                             ilojlo_rreq,    ilojlo_sreq,    ihijhi_rreq,    ihijhi_sreq,
                             ihijlo_rreq,    ihijlo_sreq,    ilojhi_rreq,    ilojhi_sreq,
                             ilojloklo_rreq, ilojloklo_sreq, ilojhiklo_rreq, ilojhiklo_sreq,
                             ihijloklo_rreq, ihijloklo_sreq, ihijhiklo_rreq, ihijhiklo_sreq,
                             ilojlokhi_rreq, ilojlokhi_sreq, ilojhikhi_rreq, ilojhikhi_sreq,
                             ihijlokhi_rreq, ihijlokhi_sreq, ihijhikhi_rreq, ihijhikhi_sreq,
                             iloklo_halo_rreq,iloklo_halo_sreq,ilokhi_halo_rreq,ilokhi_halo_sreq,
                             ihiklo_halo_rreq,ihiklo_halo_sreq,ihikhi_halo_rreq,ihikhi_halo_sreq,
                             jloklo_halo_rreq,jloklo_halo_sreq,jlokhi_halo_rreq,jlokhi_halo_sreq,
                             jhiklo_halo_rreq,jhiklo_halo_sreq,jhikhi_halo_rreq,jhikhi_halo_sreq,
                            ])
    else
        stats = MPI.Waitall!([ihi_rreq, ihi_sreq, ilo_rreq, ilo_sreq,
                              jhi_rreq, jhi_sreq, jlo_rreq, jlo_sreq,
                              khi_rreq, khi_sreq, klo_rreq, klo_sreq,
                              ])
    end

    return nothing
end

