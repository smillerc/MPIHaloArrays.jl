<h1> <img src="docs/src/assets/logo.png" alt="MPIHaloArrays.jl" width="50"> MPIHaloArrays.jl </h1>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)


MPIHaloArrays is a high-level array type to facilitate halo, or ghost-cell exchanges commonly found in large-scale PDE codes. Inspiration was taken from [`MPIArrays.jl`](https://github.com/barche/MPIArrays.jl) and [`ImplicitGlobalGrid.jl`](https://github.com/eth-cscs/ImplicitGlobalGrid.jl). Domains can be decomposed into 1, 2, or 3D parallel topologies. 

## Installation

The package can be installed with (soon to be registered)

```julia
pkg> add https://github.com/smillerc/MPIHaloArrays.jl
```

## Documentation

- [[**STABLE**](https://smillerc.github.io/MPIHaloArrays.jl/stable)] &mdash; **most recently tagged version of the documentation.**
- [[**DEV**](https://smillerc.github.io/MPIHaloArrays.jl/dev)] &mdash; **most recent development version of the documentation.**

## Basic Usage

Halo exchange is a common practice in large-scale PDE codes that decompose the domain into many sub-domains. Neighbor information is exchanged at regular intervals through "ghost" or "halo" cell regions. The image below shows an example from a 1D array that has a halo region of 3 cells.

<img src="docs/src/assets/1d_halo.png" alt="MPIHaloArrays.jl" width="600">

Halo exchanges can be done in multiple dimensions. At the moment, `MPIHaloArrays.jl` limits this to 1-3D arrays, but this will be extended in the future. The example below shows how to set up the initial array, fill halo/domain cells, do a halo exchange, and more.

**Currently arrays are limited to 1, 2, or 3D. This will be addressed in future versions**
```julia
using MPI, MPIHaloArrays

MPI.Init()
rank = MPI.Comm_rank(comm)

# Create the MPI topology
topo = CartesianTopology([4,4], # use a 4x4 decomposition
                         [true, true]) # periodic in both dimensions   

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
ilo, ihi, jlo, jhi = localindices(x) # -> useful for looping without going into halo regions

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

## Examples

- [2D Heat Diffusion](docs/examples/04-diffusion2d.jl)


[docs-stable-url]: https://smillerc.github.io/MPIHaloArrays.jl/stable
[docs-dev-url]: https://smillerc.github.io/MPIHaloArrays.jl/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg