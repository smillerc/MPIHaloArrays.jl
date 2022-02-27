
"""
Partition the array `A` on the rank `root` into chunks based on the given parallel toplogy. The array data in `A`
does not have halo regions. The MPIHaloArray constructor adds the halo regions.
This returns a `MPIHaloArray`

# Arguments
 - `A`: Global array to be split up into chunks and sent to all ranks. This does **not** include halo cells
 - `root`: MPI rank that `A` lives on
 - `nhalo`: Number of halo cells to create
"""
function scatterglobal(A::AbstractArray{T, 1}, root::Int, nhalo::Int, topology::ParallelTopology; do_corners = true, com_model = :p2p) where {T}

    if A isa Base.UnitRange A = collect(A) end

    if topology.rank == root
        # Divide the domain into sub-domains or tiles
        sizes = zeros(Int64, topology.nprocs) 
        for p in 1:topology.nprocs
            ilo, ihi = tile_indices_1d(length(A), topology.nprocs, p)
            sizes[p] = ihi-ilo+1
        end

        size_ubuf = UBuffer(sizes, 1)
        A_vbuf = VBuffer(A, sizes) 
    else
        # these variables can be set to `nothing` on non-root processes
        size_ubuf = UBuffer(nothing)
        A_vbuf = VBuffer(nothing)
    end
        
    local_size = MPI.Scatter(size_ubuf, NTuple{1, Int64}, root, topology.comm)

    A_local = MPI.Scatterv!(A_vbuf, Array{T, 1}(undef, local_size), root, topology.comm)

    return MPIHaloArray(A_local, topology, nhalo; do_corners = do_corners, com_model = com_model)
end

function scatterglobal(A::AbstractArray{T, 2}, root::Int, nhalo::Int, topology::ParallelTopology; do_corners = true, com_model = :p2p) where {T}

    if A isa Base.ReshapedArray A = collect(A) end

    if topology.rank == root
        # Divide the domain into sub-domains or tiles
        sizes = zeros(Int64, 2, topology.nprocs) 
        for p in 1:topology.nprocs
            ilo, ihi, jlo, jhi = tile_indices_2d(size(A), topology.nprocs, p)
            sizes[1,p] = ihi-ilo+1
            sizes[2,p] = jhi-jlo+1
        end

        size_ubuf = UBuffer(sizes, 2)
        counts = vec(prod(sizes, dims=1))
        A_vbuf = VBuffer(A, counts) 
    else
        # these variables can be set to `nothing` on non-root processes
        size_ubuf = UBuffer(nothing)
        A_vbuf = VBuffer(nothing)
    end
        
    local_size = MPI.Scatter(size_ubuf, NTuple{2, Int64}, root, topology.comm)

    A_local = MPI.Scatterv!(A_vbuf, Array{T, 2}(undef, local_size), root, topology.comm)

    return MPIHaloArray(A_local, topology, nhalo; do_corners = do_corners, com_model = com_model)
end

function scatterglobal(A::AbstractArray{T, 3}, root::Int, nhalo::Int, topology::ParallelTopology; do_corners = true, com_model = :p2p) where {T}

    if A isa Base.ReshapedArray A = collect(A) end

    if topology.rank == root
        # Divide the domain into sub-domains or tiles
        sizes = zeros(Int64, 3, topology.nprocs) 
        for p in 1:topology.nprocs
            ilo, ihi, jlo, jhi, klo, khi = tile_indices_3d(size(A), topology.nprocs, p)
            sizes[1,p] = ihi-ilo+1
            sizes[2,p] = jhi-jlo+1
            sizes[3,p] = khi-klo+1
        end

        size_ubuf = UBuffer(sizes, 3)
        counts = vec(prod(sizes, dims=1))
        A_vbuf = VBuffer(A, counts) 
    else
        # these variables can be set to `nothing` on non-root processes
        size_ubuf = UBuffer(nothing)
        A_vbuf = VBuffer(nothing)
    end
        
    local_size = MPI.Scatter(size_ubuf, NTuple{3, Int64}, root, topology.comm)
    A_local = MPI.Scatterv!(A_vbuf, Array{T, 3}(undef, local_size), root, topology.comm)

    return MPIHaloArray(A_local, topology, nhalo; do_corners = do_corners, com_model = com_model)
end

"""
Gather all `MPIHaloArray`s onto the `root` MPI rank and stitch together. This will ignore halo region data and create a `Array` 
that represents the global state.

# Arguments
 - `A`: MPIHaloArray
 - `root`: MPI rank to gather `A` to
"""
function gatherglobal(A::MPIHaloArray{T, 1}; root=0) where {T}

    # Index ranges excluding the halo regions
    ilo, ihi = A.local_indices[1].domain
    local_data = @view A.data[ilo:ihi]
    local_size = ihi-ilo+1

    sizes = MPI.Gather(local_size, root, A.topology.comm)

    if A.topology.rank == root
        size_ubuf = UBuffer(sizes, 1)
        output_dim = sum(sizes)
        output_vbuf = VBuffer(Array{T, 1}(undef, output_dim), sizes)
    else
        # these variables can be set to `nothing` on non-root processes
        size_ubuf = UBuffer(nothing)
        output_vbuf = VBuffer(nothing)
    end       

    A_global = MPI.Gatherv!(local_data, output_vbuf, root, A.topology.comm)

    return A_global
end

function gatherglobal(A::MPIHaloArray{T, 2}; root=0) where {T}

    # Index ranges excluding the halo regions
    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain
    local_data = @view A.data[ilo:ihi,jlo:jhi]
    local_size = (ihi-ilo+1, jhi-jlo+1)

    sizes = MPI.Gather(local_size, root, A.topology.comm)

    if A.topology.rank == root
        size_ubuf = UBuffer(sizes, 2)
        counts = [prod(s) for s in sizes]
        global_dims = A.topology.global_dims |> Tuple

        output_dim = globalsize(sizes, global_dims)
        output_vbuf = VBuffer(Array{T, 2}(undef, output_dim), counts)
    else
        # these variables can be set to `nothing` on non-root processes
        size_ubuf = UBuffer(nothing)
        output_vbuf = VBuffer(nothing)
    end       

    A_global = MPI.Gatherv!(local_data, output_vbuf, root, A.topology.comm)

    return A_global
end

function gatherglobal(A::MPIHaloArray{T, 3}; root=0) where {T}

    # Index ranges excluding the halo regions
    ilo, ihi = A.local_indices[1].domain
    jlo, jhi = A.local_indices[2].domain
    klo, khi = A.local_indices[3].domain
    local_data = @view A.data[ilo:ihi,jlo:jhi, klo:khi]
    local_size = (ihi-ilo+1, jhi-jlo+1, khi-klo+1)

    sizes = MPI.Gather(local_size, root, A.topology.comm)

    if A.topology.rank == root
        size_ubuf = UBuffer(sizes, 2)
        counts = [prod(s) for s in sizes]
        global_dims = A.topology.global_dims |> Tuple

        output_dim = globalsize(sizes, global_dims)
        output_vbuf = VBuffer(Array{T, 3}(undef, output_dim), counts)
    else
        # these variables can be set to `nothing` on non-root processes
        size_ubuf = UBuffer(nothing)
        output_vbuf = VBuffer(nothing)
    end       

    A_global = MPI.Gatherv!(local_data, output_vbuf, root, A.topology.comm)

    return A_global
end

"""Find the global dims of based on the list of local `MPIHaloArray` sizes"""
function globalsize(sizes::Vector{Tuple{T, T, T}}, globaldims::NTuple{3,T}) where {T<:Integer}
	dim_i, dim_j, dim_k = get_dims(sizes, globaldims)
	
	global_dim = (sum(dim_i[1:globaldims[1],1,1]), 
                  sum(dim_j[1,1:globaldims[2],1]), 
                  sum(dim_k[1,1,1:globaldims[3]]))
	
    return global_dim
end

function globalsize(sizes::Vector{Tuple{T, T}}, globaldims::NTuple{2,T}) where {T<:Integer}
	dim_i, dim_j = get_dims(sizes, globaldims)
	
	global_dim = (sum(dim_i[1:globaldims[1],1]), sum(dim_j[1,1:globaldims[2]]))
    return global_dim
end

"""Get the dimensions of each chunk"""
function get_dims(sizes::Vector{Tuple{T, T, T}}, globaldims::NTuple{3,T}) where {T<:Integer}

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

function get_dims(sizes::Vector{Tuple{T, T}}, globaldims::NTuple{2,T}) where {T<:Integer}

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