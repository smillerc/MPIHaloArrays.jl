
"""
Helper functions to get the low side halo and domain starting/ending indices

# Arguments
 - `field`: Array 
 - `dim`: Dimension to check indices on
 - `nhalo`: Number of halo entries
"""
@inline function lo_indices(field, dim, nhalo)
    CI = CartesianIndices(field)
    lo_indices = first(CI) |> Tuple
    hi_indices = last(CI) |> Tuple
    
    lo_halo_start = lo_indices[dim]        # start index of the halo region

    lo_dom_start = nhalo + 1 # start index of the inner domain
    lo_dom_end = 2 * nhalo   # end index of the inner domain
    lo_halo_end = lo_halo_start + nhalo - 1 # end index of the halo region

    return (lo_halo_start, lo_halo_end, lo_dom_start, lo_dom_end)
end

"""Helper functions to get the high side halo and domain starting/ending indices"""
@inline function hi_indices(field, dim, nhalo)

    if dim > length(size(field)) error("dim > length(size(field))") end
    hi_halo_end = size(field, dim)       # end index of the halo region
    hi_halo_start = hi_halo_end - nhalo + 1  # start index of the halo region

    hi_dom_end = hi_halo_start - 1   # end index of the inner domain
    hi_dom_start = hi_dom_end - nhalo + 1 # start index of the inner domain

    return (hi_dom_start, hi_dom_end, hi_halo_start, hi_halo_end)
end