# MPIHaloArrays

MPIHaloArrays is a high-level array type to help with halo, or ghost-cell exchanges commonly found in large-scale PDE problems. Very similar in goals and design to `MPIArrays.jl` and `ImplicitGlobalGrid.jl`. Domains can be decomposed into 1, 2, or 3D topology. Currently arrays are limited to 1, 2, or 3D.

## Installation

The package is registered and can be installed with

```julia
pkg> add MPIHaloArrays
```

## Basic Usage

This example shows how to set up the initial array, fill the halo/domain cells, and do a halo exchange
```julia
using MPI, MPIHaloArrays

MPI.Init()
rank = MPI.Comm_rank(comm)

# Create the MPI topology
topo = CartesianTopology([4,4], # use a 4x4 decomposition
                         [true, true]; # periodic in both dimensions   
                         do_corners=false) # exchange corner halo regions (significant speed advantage if you don't need it)

nhalo = 2 # Number of halo cells in each dimension (fixed for all dimensions)
N = 200

# create the array type; this pads the data on all sides with halo regions
x = MPIHaloArray(rand(N,N), topo, nhalo)

# fill all the halo regions with -1
fillhalo!(A, -1)

# fill the domain region with the current rank
filldomain!(A, rank)

# local (current rank) indexing works just like a normal array
x[1,1] .= 2.0

# Get the local/global indices of the _domain_ data (not including the halo cells)
ilo, ihi, jlo, jhi = local_domain_indices(x) # -> useful for looping without going into halo regions
ilo_g, ihi_g, jlo_g, jhi_g = global_domain_indices(x)

# Exchange data with neighbors
updatehalo!(x)

GC.gc()
MPI.Finalize()
```

Scatter and gather operations are also defined with `scatterglobal` and `gatherglobal`.

```julia
rank = 0 # MPI rank to scatter from / gather to

# start with a global Base.Array type to decompose and scatter to each rank
ni = 512; nj = 256
A_global = reshape(1:ni*nj, ni, nj);

# scatter - this internally converts A_global to multiple halo arrays. This is why
#           the nhalo and topology types are needed
A_local = scatterglobal(A_global, root, nhalo, topology) # -> returns a MPIHaloArray

# do some work...

# and now gather the decomposed domain and store on the root rank of choice
A_global_result = gatherglobal(A_local; root=root) # -> returns a Base.Array
```


At the moment, reductions are not implemented, but will be in the future...