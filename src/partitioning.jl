"""
    split_count(N::Integer, n::Integer)

Return a vector of `n` integers which are approximately equally sized and sum to `N`.
"""
function split_count(N::Integer, n::Integer)
    q,r = divrem(N, n)
    return [i <= r ? q+1 : q for i = 1:n]
end

"""
	get_subdomain_dimension_sizes(A, tile_dims, A_halo_dims)

Get the size along each dimension in `(i,j,k)` of the subdomain, based on a given array `A`. 
The `tile_dims` is the shape of the global domain, e.g., (4,2) means 4 tiles or subdomains in 
`i` and 2 in `j`. `A_halo_dims` is the tuple of which dimensions the halo exchanges take place on, e.g. `(2,3)`.

# Example
```julia
A = rand(4,200,100); dims=(2,3), tile_dims=(4,2)
get_subdomain_dimension_sizes(A, tile_dims, dims) # [[i][j]] --> [[50,50,50,50],[100,100]]
```
"""
function get_subdomain_dimension_sizes(A_size, tile_dims, A_halo_dims)
	tile_sizes = Vector{Vector{Int}}(undef, 0)
	new_tile_dims = match_tile_halo_dim_sizes(tile_dims, A_halo_dims)
	for (nd, nt) in zip(A_halo_dims, new_tile_dims)
		sz = split_count(A_size[nd], nt)
		push!(tile_sizes,sz)
	end
	tile_sizes
end

"""
	match_tile_halo_dim_sizes(tile_dims, halo_dims)

Ensure that the tile dimension tuple is the same length as the halo dim tuple. If not, then pad with ones.
"""
function match_tile_halo_dim_sizes(tile_dims, halo_dims)
	new_tile_dims = ones(Int,length(halo_dims))
	for i in 1:length(tile_dims)
		new_tile_dims[i] = tile_dims[i]
	end
	new_tile_dims |> Tuple
end

"""
Get the size of each subdomain given the tile dimensions and number of halo cells
"""
function get_subdomain_sizes(A_size, tile_dims::NTuple{1,Int}, halo_dims)

	# the array that hold sizes of each subdomain; [size, proc_id]
	sizes = zeros(Int, length(halo_dims), prod(tile_dims))
	
	sub_sizes = get_subdomain_dimension_sizes(A_size, tile_dims, halo_dims)
	ranges = [UnitRange(1,length(s)) for s in sub_sizes] |> Tuple

	for i in LinearIndices(ranges)
		sizes[1,i] = sub_sizes[1][i]
	end

	sizes	
end

function get_subdomain_sizes(A_size, tile_dims::NTuple{2,Int}, halo_dims)

	# the array that hold sizes of each subdomain; [size, proc_id]
	sizes = zeros(Int, length(halo_dims), prod(tile_dims))
	
	sub_sizes = get_subdomain_dimension_sizes(A_size, tile_dims, halo_dims)
	ranges = [UnitRange(1,length(s)) for s in sub_sizes] |> Tuple

	LI = LinearIndices(ranges)
	for Idx in CartesianIndices(LI)
		i, j = Tuple(Idx)
		I = LI[Idx]
		
		sizes[1,I] = sub_sizes[1][i]
		sizes[2,I] = sub_sizes[2][j]
	end
	sizes	
end

function get_subdomain_sizes(A_size, tile_dims::NTuple{3,Int}, halo_dims)

	# the array that hold sizes of each subdomain; [[dimi, dimj, etc.], proc_id]
	sizes = zeros(Int, length(halo_dims), prod(tile_dims))
	
	sub_sizes = get_subdomain_dimension_sizes(A_size, tile_dims, halo_dims)
	ranges = [UnitRange(1,length(s)) for s in sub_sizes] |> Tuple

	LI = LinearIndices(ranges)
	for Idx in CartesianIndices(LI)
		i, j, k = Tuple(Idx)
		I = LI[Idx]
		
		sizes[1,I] = sub_sizes[1][i]
		sizes[2,I] = sub_sizes[2][j]
		sizes[3,I] = sub_sizes[3][k]
	end
	sizes	
end

"""
Get the starting and ending indices based on the size of each subdomain. These are
the global ilo/ihi values.
"""
function get_istarts_ends(isizes)
	istarts = similar(isizes)
	iends = similar(isizes)
	for i in 1:length(isizes)
		if i > 1
			istarts[i] = istarts[i-1] + isizes[i-1]
			iends[i] = istarts[i] + isizes[i]-1
		else
			istarts[i] = 1
			iends[i] = isizes[i]
		end
	end
	return istarts, iends
end

"""
Get the global indices for each subdomain. This is the tuple of lo and hi indices for each dimension
"""
function get_subdomain_indices(A_size, domain_size::NTuple{3,Int}, halo_dims)

	ndims = length(domain_size)
	isizes, jsizes, ksizes = get_subdomain_dimension_sizes(A_size, domain_size, halo_dims)
	iloihi = get_istarts_ends(isizes)
	jlojhi = get_istarts_ends(jsizes)
	klokhi = get_istarts_ends(ksizes)
	lohi_indices = Vector{NTuple{ndims*2, Int}}(undef, 0)
	for k in 1:domain_size[3]
		for j in 1:domain_size[2]
			for i in 1:domain_size[1]
				ilo = iloihi[1][i]
				ihi = iloihi[2][i]
				jlo = jlojhi[1][j]
				jhi = jlojhi[2][j]
				klo = klokhi[1][k]
				khi = klokhi[2][k]
				push!(lohi_indices, (ilo,ihi,jlo,jhi,klo,khi)) 
			end
		end
	end
	lohi_indices
end

function get_subdomain_indices(A_size, domain_size::NTuple{2,Int}, halo_dims)

	ndims = length(domain_size)
	isizes, jsizes = get_subdomain_dimension_sizes(A_size, domain_size, halo_dims)
	iloihi = get_istarts_ends(isizes)
	jlojhi = get_istarts_ends(jsizes)
	lohi_indices = Vector{NTuple{ndims*2, Int}}(undef, 0)
	for j in 1:domain_size[2]
			for i in 1:domain_size[1]
			ilo = iloihi[1][i]
			ihi = iloihi[2][i]
			jlo = jlojhi[1][j]
			jhi = jlojhi[2][j]
			push!(lohi_indices, (ilo,ihi,jlo,jhi)) 
		end
	end
	lohi_indices
end

function get_subdomain_indices(A_size, domain_size::NTuple{1,Int}, halo_dims)

	ndims = length(domain_size)
	isizes = get_subdomain_dimension_sizes(A_size, domain_size, halo_dims)
	iloihi = get_istarts_ends(isizes[1])
	lohi_indices = Vector{NTuple{ndims*2, Int}}(undef, 0)
	for i in 1:domain_size[1]
		ilo = iloihi[1][i]
		ihi = iloihi[2][i]
		push!(lohi_indices, (ilo,ihi)) 
	end
	lohi_indices
end