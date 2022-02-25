using LinearAlgebra: norm

"""Return all common denominators of n"""
function denominators(n::Integer)
    denominators = Vector{Int}(undef, 0)
    for i in 1:n
        if mod(n, i) == 0
            push!(denominators, i)
        end
    end
    return denominators
end

"""Returns the optimal number of tiles in (i,j) given total number of tiles n"""
function num_2d_tiles(n)
    # find all common denominators of the total number of images
    denoms = denominators(n)

    # find all combinations of common denominators
    # whose product equals the total number of images
    dim1 = Vector{Int}(undef, 0)
    dim2 = Vector{Int}(undef, 0)
    for j in 1:length(denoms)
        for i in 1:length(denoms)
            if denoms[i] * denoms[j] == n
                push!(dim1, denoms[i])
                push!(dim2, denoms[j])
            end
        end
    end
    # pick the set of common denominators with the minimal norm
    # between two elements -- rectangle closest to a square
    num_2d_tiles = [dim1[1], dim2[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i]] .- sqrt(n))
        n2 = norm(num_2d_tiles .- sqrt(n))
        if n1 < n2
            num_2d_tiles = [dim1[i], dim2[i]]
        end
    end
    return num_2d_tiles
end

"""Returns the optimal number of tiles in (i,j,k) given total number of tiles n"""
function num_3d_tiles(n)
    # find all common denominators of the total number of images
    denoms = denominators(n)

    # find all combinations of common denominators
    # whose product equals the total number of images
    dim1 = Vector{Int}(undef, 0)
    dim2 = Vector{Int}(undef, 0)
    dim3 = Vector{Int}(undef, 0)
    for k in 1:length(denoms)
        for j in 1:length(denoms)
            for i in 1:length(denoms)
                if denoms[i] * denoms[j] * denoms[k] == n
                    push!(dim1, denoms[i])
                    push!(dim2, denoms[j])
                    push!(dim3, denoms[k])
                end
            end
        end
    end
    # pick the set of common denominators with the minimal norm
    # between two elements -- rectangle closest to a square
    num_3d_tiles = [dim1[1], dim2[1], dim3[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i], dim3[i]] .- sqrt(n))
        n2 = norm(num_3d_tiles .- sqrt(n))
        if n1 < n2
            num_3d_tiles = [dim1[i], dim2[i], dim3[i]]
        end
    end
    return num_3d_tiles
end

"""
Given an input I dimensions of the total computational domain,
returns an array of start and end indices [ilo,ihi]
"""
function tile_indices_1d(dims::Integer, ntiles::Integer, id::Integer)
    indices = zeros(Int, 2)
    tile_size = dims ÷ ntiles

    # start and end indices assuming equal tile sizes
    indices[1] = (id - 1) * tile_size + 1
    indices[2] = indices[1] + tile_size - 1

    # if we have any remainder, distribute it to the tiles at the end
    offset = ntiles - mod(dims, ntiles)
    if id > offset
        indices[1] = indices[1] + id - offset - 1
        indices[2] = indices[2] + id - offset
    end
    return indices
end

"""
Given an input (I,J) dimensions of the total computational domain,
returns an array of start and end indices [ilo,ihi,jlo,jhi]
"""
function tile_indices_2d(dims, ntiles::Integer, id::Integer)
    indices = zeros(Int, 4)
    tiles = num_2d_tiles(ntiles)
    tiles_ij = tile_id_to_ij(id, ntiles)
    indices[1:2] = tile_indices_1d(dims[1], tiles[1], tiles_ij[1])
    indices[3:4] = tile_indices_1d(dims[2], tiles[2], tiles_ij[2])
    return indices
end

"""
Given an input (I,J,K) dimensions of the total computational domain,
returns an array of start and end indices [ilo,ihi,jlo,jhi,klo,khi]
"""
function tile_indices_3d(dims, ntiles::Integer, id::Integer)
    indices = zeros(Int, 6)
    tiles = num_3d_tiles(ntiles)
    tiles_ij = tile_id_to_ijk(id, ntiles)
    indices[1:2] = tile_indices_1d(dims[1], tiles[1], tiles_ij[1])
    indices[3:4] = tile_indices_1d(dims[2], tiles[2], tiles_ij[2])
    indices[5:6] = tile_indices_1d(dims[3], tiles[3], tiles_ij[3])
    return indices
end

"""Given tile id in a 1D layout, returns the corresponding tile indices in a 2D layout"""
function tile_id_to_ij(id::Integer, ntiles::Integer)
    if id < 1
        @error("Invalid tile id")
    end
    I, J = num_2d_tiles(ntiles)
    CI = CartesianIndices((1:I, 1:J))
    ij = Tuple(CI[id])
    return ij
end

"""Given tile id in a 1D layout, returns the corresponding tile indices in a 3D layout"""
function tile_id_to_ijk(id::Integer, ntiles::Integer)
    if id < 1
        @error("Invalid tile id")
    end
    I, J, K = num_3d_tiles(ntiles)
    CI = CartesianIndices((1:I, 1:J, 1:K))
    ijk = Tuple(CI[id])
    return ijk
end