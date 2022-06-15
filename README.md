<h1> <img src="docs/src/assets/logo.png" alt="MPIHaloArrays.jl" width="50"> MPIHaloArrays.jl </h1>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/MPIHaloArrays)](https://pkgs.genieframework.com?packages=MPIHaloArrays).




MPIHaloArrays provides a high-level array type to facilitate halo, or ghost-cell exchanges commonly found in large-scale PDE codes. The `MPIHaloArray` type is a subtype of `AbstractArray`.

Inspiration was taken from [`MPIArrays.jl`](https://github.com/barche/MPIArrays.jl) and [`ImplicitGlobalGrid.jl`](https://github.com/eth-cscs/ImplicitGlobalGrid.jl). Domains can be decomposed into 1, 2, or 3D parallel topologies. 

## Installation

The package can be installed with

```julia
pkg> add MPIHaloArrays
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
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)
const root = 0

# Create the MPI topology
topo = CartesianTopology(comm, [4,4], # use a 4x4 decomposition
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
A[1,1] .= 2.0

# Get the local/global indices of the _domain_ data (not including the halo cells)
ilo, ihi, jlo, jhi = local_domain_indices(x) # -> useful for looping without going into halo regions

# Exchange data with neighbors
updatehalo!(x)

GC.gc()
MPI.Finalize()
```

Scatter and gather operations are also defined with `scatterglobal` and `gatherglobal`.

```julia
root = 0 # MPI rank to scatter from / gather to

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

## Setting the halo exchange dimensions

In some cases, you may want a multi-dimensional array that only does the halo exchange on a subset of the dimensions. For example, an array `U`, has dimensions `[[ρ, u, v, p], i, j]`, where `[i,j]` represent 2D grid coordinates and `[ρ, u, v, p]` are the density, x/y velocity, and pressure at each `[i,j]` coordinate. This can be done with the following syntax:

```julia
using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)
const root = 0

U = zeros(4,128,128)

# Only 1 MPI rank, but the topology is 2D
topo = CartesianTopology(comm, (1,1), (true,true))

nhalo = 2
halo_dims = (2,3)
A = MPIHaloArray(U, topo, nhalo, halo_dims)

# This will fail b/c U is 3D and the topology is only 2D; you must
# specify the halo exchange in 2D
A = MPIHaloArray(U, topo, nhalo) # -> ERROR: Mismatched topology dimensionality (2D) and halo region dimensions (3D)


```

**Note:** The default behavior selects __all__ dimensions to exchange halo data. You must provide the `halo_dims` tuple to override this.


## Interoperability

Add physical units via `Unitful.jl`
```julia
using MPIHaloArrays, Unitful
data = rand(10,10) * u"m"
A = MPIHaloArray(data, topology, 2)
```

Add uncertainty via `Measurements.jl`
```julia
using MPIHaloArrays, Unitful, Measurements
data = (rand(10,10) .± 0.1) * u"m"
A = MPIHaloArray(data, topology, 2)
```


## Examples

A slightly more useful example that performs 2D heat diffusion is shown [here](docs/examples/04-diffusion2d.jl). This shows how to
 - Scatter initial conditions from the root node to each MPI process with `scatterglobal()`
 - Perform a stencil operation within the current `MPIHaloArray`. This looks like any other normal array loop, but the bounds are determined by the `MPIHaloArray` via `local_domain_indices()`
 - Update halo cells / neighbor information. Periodic boundary conditions are also handled by the `CartesianTopology` type.
 - Gather results to the root node for plotting/output with `gatherglobal()`


## Exported functions/types

- `MPIHaloArray`: An array type that extends `AbstractArray` to provide MPI neighbor communication for halo or ghost cells
- `AbstractParallelTopology`, `CartesianTopology`: MPI Topology types to manage neighbor information
- `neighbor(), neighbors()`, `[i,j,k]lo_neighbor()`, `[i,j,k]hi_neighbor()`: Extract neighbors of the current MPI rank
- `lo_indices()`, `hi_indices()`: Local indices of the current MPIHaloArray. Used for loop limits that ignore halo regions
- `fillhalo!()`: Fill the halo cells with a scalar value
- `filldomain!()`: Fill the domain cells with a scalar value
- `updatehalo!()`: Perform neighbor communication / halo exchange
- `scatterglobal()`: Distribute/scatter a global array to multiple ranks - returns a local `MPIHaloArray` for each rank
- `gatherglobal()`: Gather `MPIHaloArray`s to a root MPI rank - returns an `AbstractArray` on the root node

[docs-stable-url]: https://smillerc.github.io/MPIHaloArrays.jl/stable
[docs-dev-url]: https://smillerc.github.io/MPIHaloArrays.jl/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
