# MPIHaloArrays.jl - Work in progress


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/dev)

A high-level array type to help with halo, or ghost-cell exchanges commonly found in large-scale PDE problems. Very similar in goals and design to `MPIArrays.jl` and `ImplicitGlobalGrid.jl`.

```julia
using MPI, MPIHaloArrays

MPI.Init()

# Create the MPI topology, which here is a Cartesian 4x4 domain (using 16 cores)
topo = CartesianTopology([4,4], [true, true])
rank = MPI.Comm_rank(MPI.COMM_WORLD)

# How many halo cells in each dimension (fixed for all dimensions)
nhalo = 2
N = 200

local_data = rand(N,N) # does not include halo cells
x = MPIHaloArray{Float64}(local_data,topo,nhalo)

updatehalo!(x)

GC.gc()
MPI.Finalize()
```