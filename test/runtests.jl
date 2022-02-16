using Test, MPI

nprocs_str = get(ENV, "JULIA_MPI_TEST_NPROCS","")

include("unit/run_unittests.jl")