
"""
Partition the array `A` on the rank `root` into chunks based on the given parallel toplogy. The array data in `A`
does not have halo regions. The MPIHaloArray constructor adds the halo regions.
This returns a `MPIHaloArray`

# Arguments
 - `A`: Global array to be split up into chunks and sent to all ranks. This does **not** include halo cells
 - `root`: MPI rank that `A` lives on
 - `nhalo`: Number of halo cells to create
 - `halo_dims`: Tuple of the dimensions that halo exchanges occur on (not fully working yet)
"""
function scatterglobal(A::AbstractArray{T, 1}, root::Int, nhalo::Int,
    topology::AbstractParallelTopology; halo_dims = (1), do_corners = true, com_model = :p2p) where {T, N}

    if A isa Base.UnitRange A = collect(A) end

    # coords = MPI.Allgather(topology.coords, topology.comm) # -> defaults to a NTuple{3,Int}
    # subdomain_coords = [c[1:topology.dimension] for c in coords] # strip the unused dimensions
    global_domain_size = size(topology)

    sizes = get_subdomain_sizes(size(A), global_domain_size, halo_dims)
    indices = get_subdomain_indices(size(A), global_domain_size, halo_dims)

    remote_size = sizes[:,topology.rank + 1] |> Tuple

    A_local = zeros(eltype(A), remote_size)
    remote_buf = MPI.Buffer(A_local)

    reqs = Vector{MPI.Request}(undef, 0)

    if topology.rank == root
        for sendrank in 0:topology.nprocs-1

            # Get the indices on the root buffer to send to the remote buffer
            ilo, ihi = indices[sendrank + 1]
            data_on_root = @view A[ilo:ihi]
            root_buf = MPI.Buffer(data_on_root)
            sendtag = sendrank + 1000
            sreq =  MPI.Isend(root_buf, sendrank, sendtag, topology.comm)
            push!(reqs, sreq)
        end
    end

    recievetag = topology.rank + 1000
    rreq = MPI.Irecv!(remote_buf, root, recievetag, topology.comm)
    push!(reqs, rreq)

    MPI.Waitall!(reqs)
    return MPIHaloArray(A_local, topology, nhalo; do_corners = do_corners, com_model = com_model)
end

function scatterglobal(A::AbstractArray{T, 2}, root::Int, nhalo::Int,
    topology::AbstractParallelTopology; halo_dims = (1, 2), do_corners = true, com_model = :p2p) where {T, N}

    if A isa Base.ReshapedArray A = collect(A) end

    # coords = MPI.Allgather(topology.coords, topology.comm) # -> defaults to a NTuple{3,Int}
    # subdomain_coords = [c[1:topology.dimension] for c in coords] # strip the unused dimensions
    global_domain_size = size(topology)

    sizes = get_subdomain_sizes(size(A), global_domain_size, halo_dims)
    indices = get_subdomain_indices(size(A), global_domain_size, halo_dims)

    remote_size = sizes[:,topology.rank + 1] |> Tuple

    A_local = zeros(eltype(A), remote_size)
    remote_buf = MPI.Buffer(A_local)

    reqs = Vector{MPI.Request}(undef, 0)

    if topology.rank == root
        for sendrank in 0:topology.nprocs-1

            # Get the indices on the root buffer to send to the remote buffer
            ilo, ihi, jlo, jhi = indices[sendrank + 1]
            data_on_root = @view A[ilo:ihi, jlo:jhi]
            root_buf = MPI.Buffer(data_on_root)
            sendtag = sendrank + 1000
            sreq =  MPI.Isend(root_buf, sendrank, sendtag, topology.comm)
            push!(reqs, sreq)
        end
    end

    recievetag = topology.rank + 1000
    rreq = MPI.Irecv!(remote_buf, root, recievetag, topology.comm)
    push!(reqs, rreq)

    MPI.Waitall!(reqs)
    return MPIHaloArray(A_local, topology, nhalo; do_corners = do_corners, com_model = com_model)
end

function scatterglobal(A::AbstractArray{T, 3}, root::Int, nhalo::Int,
    topology::AbstractParallelTopology; halo_dims = (1, 2, 3), do_corners = true, com_model = :p2p) where {T, N}

    if A isa Base.ReshapedArray A = collect(A) end

    # coords = MPI.Allgather(topology.coords, topology.comm) # -> defaults to a NTuple{3,Int}
    # subdomain_coords = [c[1:topology.dimension] for c in coords] # strip the unused dimensions
    global_domain_size = size(topology)

    sizes = get_subdomain_sizes(size(A), global_domain_size, halo_dims)
    indices = get_subdomain_indices(size(A), global_domain_size, halo_dims)

    remote_size = sizes[:,topology.rank + 1] |> Tuple

    A_local = zeros(eltype(A), remote_size)
    remote_buf = MPI.Buffer(A_local)

    reqs = Vector{MPI.Request}(undef, 0)

    if topology.rank == root
        for sendrank in 0:topology.nprocs-1

            # Get the indices on the root buffer to send to the remote buffer
            ilo, ihi, jlo, jhi, klo, khi = indices[sendrank + 1]
            # @show sendrank, (ilo, ihi, jlo, jhi, klo, khi)
            data_on_root = @view A[ilo:ihi, jlo:jhi, klo:khi]
            root_buf = MPI.Buffer(data_on_root)
            sendtag = sendrank + 1000
            sreq =  MPI.Isend(root_buf, sendrank, sendtag, topology.comm)
            push!(reqs, sreq)
        end
    end

    recievetag = topology.rank + 1000
    rreq = MPI.Irecv!(remote_buf, root, recievetag, topology.comm)
    push!(reqs, rreq)

    MPI.Waitall!(reqs)
    return MPIHaloArray(A_local, topology, nhalo; do_corners = do_corners, com_model = com_model)
end

"""
Gather all `MPIHaloArray`s onto the `root` MPI rank and stitch together. This will ignore halo region data and create a `Array`
that represents the global state.

# Arguments
 - `A`: MPIHaloArray
 - `root`: MPI rank to gather `A` to
 - `halo_dims`: Tuple of the dimensions that halo exchanges occur on (not fully working yet)
"""
function gatherglobal(A::MPIHaloArray{T, N, 1}; root=0) where {T, N}

    halo_dims = A.halo_dims

    # Index ranges excluding the halo regions
    ilo_l, ihi_l = localindices(A)
    local_indices = (ilo_l, ihi_l)
    local_data = @view A.data[ilo_l:ihi_l]
    local_size = size(local_data)
    local_buf = MPI.Buffer(local_data)

    # Put all the indices on root so the gather below can use them
    # indices = MPI.Gather(local_indices, root, A.topology.comm)
    global_domain_size = size(A.topology)
    sizes = MPI.Gather(local_size, root, A.topology.comm)

    # Allocate the data on root
    if A.topology.rank == root
        global_dims = size(A.topology)
        output_dim = sum(s[1] for s in sizes)
        A_global = Array{T}(undef, output_dim)
        indices = get_subdomain_indices(size(A_global), global_domain_size, halo_dims)
    end

    reqs = Vector{MPI.Request}(undef, 0)

    # All ranks send to root
    sendtag = A.topology.rank + 1000
    sreq =  MPI.Isend(local_buf, root, sendtag, A.topology.comm)
    push!(reqs, sreq)

    # Loop through all ranks and recieve from each
    if A.topology.rank == root
        for receiverank in 0:A.topology.nprocs-1

            receivetag = receiverank + 1000
            ilo, ihi = indices[receiverank + 1]
            data_on_root = @view A_global[ilo:ihi]
            root_buf = MPI.Buffer(data_on_root)
            rreq = MPI.Irecv!(root_buf, receiverank, receivetag, A.topology.comm)
            push!(reqs, rreq)
        end
    else
        A_global = nothing
    end

    MPI.Waitall!(reqs)
    return A_global
end

function gatherglobal(A::MPIHaloArray{T, N, 2}; root=0) where {T, N}

    halo_dims = A.halo_dims

    # Index ranges excluding the halo regions
    ilo_l, ihi_l, jlo_l, jhi_l = localindices(A)
    local_indices = (ilo_l, ihi_l, jlo_l, jhi_l)
    local_data = @view A.data[ilo_l:ihi_l, jlo_l:jhi_l]
    local_size = size(local_data)
    local_buf = MPI.Buffer(local_data)
    # println("Rank $(A.topology.rank), local_indices: $(local_indices)\n")

    # Put all the indices on root so the gather below can use them
    # indices = MPI.Gather(local_indices, root, A.topology.comm)
    global_domain_size = size(A.topology)
    sizes = MPI.Gather(local_size, root, A.topology.comm)

    # Allocate the data on root
    if A.topology.rank == root
        global_dims = size(A.topology)
        output_dim = globalsize(sizes, global_dims)
        A_global = Array{T}(undef, output_dim)
        indices = get_subdomain_indices(size(A_global), global_domain_size, halo_dims)
    end

    reqs = Vector{MPI.Request}(undef, 0)

    # All ranks send to root
    sendtag = A.topology.rank + 1000
    sreq =  MPI.Isend(local_buf, root, sendtag, A.topology.comm)
    push!(reqs, sreq)

    # Loop through all ranks and recieve from each
    if A.topology.rank == root
        for receiverank in 0:A.topology.nprocs-1

            receivetag = receiverank + 1000
            ilo, ihi, jlo, jhi = indices[receiverank + 1]
            data_on_root = @view A_global[ilo:ihi, jlo:jhi]
            root_buf = MPI.Buffer(data_on_root)
            rreq = MPI.Irecv!(root_buf, receiverank, receivetag, A.topology.comm)
            push!(reqs, rreq)
        end
    else
        A_global = nothing
    end

    MPI.Waitall!(reqs)
    return A_global
end

function gatherglobal(A::MPIHaloArray{T, N, 3}; root=0) where {T, N}

    halo_dims = A.halo_dims

    # Index ranges excluding the halo regions
    ilo_l, ihi_l, jlo_l, jhi_l, klo_l, khi_l = localindices(A)
    local_indices = (ilo_l, ihi_l, jlo_l, jhi_l, klo_l, khi_l)
    local_data = @view A.data[ilo_l:ihi_l, jlo_l:jhi_l, klo_l:khi_l]
    local_size = size(local_data)
    local_buf = MPI.Buffer(local_data)
    # println("Rank $(A.topology.rank), local_indices: $(local_indices)\n")

    # Put all the indices on root so the gather below can use them
    # indices = MPI.Gather(local_indices, root, A.topology.comm)
    global_domain_size = size(A.topology)
    sizes = MPI.Gather(local_size, root, A.topology.comm)

    # Allocate the data on root
    if A.topology.rank == root
        global_dims = size(A.topology)
        output_dim = globalsize(sizes, global_dims)
        A_global = Array{T}(undef, output_dim)
        indices = get_subdomain_indices(size(A_global), global_domain_size, halo_dims)
    end

    reqs = Vector{MPI.Request}(undef, 0)

    # All ranks send to root
    sendtag = A.topology.rank + 1000
    sreq =  MPI.Isend(local_buf, root, sendtag, A.topology.comm)
    push!(reqs, sreq)

    # Loop through all ranks and recieve from each
    if A.topology.rank == root
        for receiverank in 0:A.topology.nprocs-1

            receivetag = receiverank + 1000
            ilo, ihi, jlo, jhi, klo, khi = indices[receiverank + 1]
            data_on_root = @view A_global[ilo:ihi, jlo:jhi, klo:khi]
            root_buf = MPI.Buffer(data_on_root)
            rreq = MPI.Irecv!(root_buf, receiverank, receivetag, A.topology.comm)
            push!(reqs, rreq)
        end
    else
        A_global = nothing
    end

    MPI.Waitall!(reqs)
    return A_global
end

"""Find the global dims of based on the list of local `MPIHaloArray` sizes"""
function globalsize(sizes::Vector{Tuple{T, T, T}}, globaldims) where {T<:Integer}
	dim_i, dim_j, dim_k = get_dims(sizes, globaldims)

	global_dim = (sum(dim_i[1:globaldims[1],1,1]),
                  sum(dim_j[1,1:globaldims[2],1]),
                  sum(dim_k[1,1,1:globaldims[3]]))

    return global_dim
end

function globalsize(sizes::Vector{Tuple{T, T}}, globaldims) where {T<:Integer}
	dim_i, dim_j = get_dims(sizes, globaldims)

	global_dim = (sum(dim_i[1:globaldims[1],1]), sum(dim_j[1,1:globaldims[2]]))
    return global_dim
end

"""Get the dimensions of each chunk"""
function get_dims(sizes::Vector{Tuple{T, T, T}}, globaldims) where {T<:Integer}

	di = reshape([s[1] for s in sizes], globaldims)
    dj = reshape([s[2] for s in sizes], globaldims)
	dk = reshape([s[3] for s in sizes], globaldims)

	for i in 1:size(di,1)
		if length(unique(di[i,:,:])) != 1
			println("=======================")
			@error("Mismatched i dimensions")
			display(di)
			println("\n=======================")
			break
		end
	end

	for j in 1:size(dj,2)
		if length(unique(dj[:,j,:])) != 1
		    println("=======================")
			@error("Mismatched j dimensions")
			display(dj)
			println("\n=======================")
			break
		end
	end

	for k in 1:size(dk,3)
		if length(unique(dk[:,:,k])) != 1
		    println("=======================")
			@error("Mismatched k dimensions")
			display(dk)
			println("\n=======================")
			break
		end
	end
	return di, dj, dk
end

function get_dims(sizes::Vector{Tuple{T, T}}, globaldims) where {T<:Integer}

	di = reshape([s[1] for s in sizes], globaldims)
    dj = reshape([s[2] for s in sizes], globaldims)

	for i in 1:size(di,1)
		if length(unique(di[i,:])) != 1
			println("=======================")
			@error("Mismatched i dimensions")
			display(di)
			println("\n=======================")
			break
		end
	end

	for j in 1:size(dj,2)
		if length(unique(dj[:,j])) != 1
		    println("=======================")
			@error("Mismatched j dimensions")
			display(dj)
			println("\n=======================")
			break
		end
	end
	return di, dj
end
