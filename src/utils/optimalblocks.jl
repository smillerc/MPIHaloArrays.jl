
using LinearAlgebra

"""Returns all common denominators of n"""
function denominators(n)
    denoms = Vector{Int}(undef, 0)
    for i in 1:n
        if mod(n, i) == 0
            push!(denoms, i)
        end
    end
    return denoms
end

"""
Returns the optimal number of blocks in 2 dimensions given total number of blocks n.

# Example

```julia
num_2d_blocks(6) = [3, 2]
```
"""
function num_2d_blocks(n)

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
    n2d_blocks = [dim1[1], dim2[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i]] .- sqrt(n))
        n2 = norm(n2d_blocks .- sqrt(n))
        if n1 < n2
            n2d_blocks[1] = dim1[i]
            n2d_blocks[2] = dim2[i]
        end
    end

    return n2d_blocks
end

function num_3d_blocks(n)
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
    n3d_blocks = [dim1[1], dim2[1], dim3[1]]
    for i in 3:length(dim1)
        n1 = norm([dim1[i], dim2[i], dim3[i]] .- sqrt(n))
        n2 = norm(n3d_blocks .- sqrt(n))
        if n1 < n2
            n3d_blocks[1] = dim1[i]
            n3d_blocks[2] = dim2[i]
            n3d_blocks[3] = dim3[i]
        end
    end

    return n3d_blocks
end
