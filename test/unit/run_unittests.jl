using Test, MPI, BenchmarkTools

testdir = @__DIR__

# No MPI needed for these tests
@time @testset "Partitioning Tests" begin include("test_partitioning.jl") end

mpi_testfiles = [
    # ("test_halo_1d.jl", 8),
    # ("test_halo_2d.jl", 8),
    # ("test_halo_3d.jl", 8),
    # ("test_indexing.jl", 16),
    # ("test_topology_1d.jl", 16),
    # ("test_topology_2d.jl", 16),
    # ("test_topology_3d.jl", 16),
    # ("test_mpihaloarray.jl", 8),
    ("test_scattergather.jl", 8),
    ]

# Run each test with the mpiexec command -- this requires special treatment due to MPI
@testset "$(f[1])" for f in mpi_testfiles
    file, nprocs = f
    mpiexec() do cmd
        # e.g. run `mpiexec -n 4 julia test_somthing.jl`
        run(`$cmd -n $nprocs $(Base.julia_cmd()) $(joinpath(testdir, file))`)
        @test true
    end
end

nothing