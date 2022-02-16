using Test, MPI

nprocs_str = get(ENV, "JULIA_MPI_TEST_NPROCS","")


testdir = @__DIR__
istest(f) = endswith(f, ".jl") && startswith(f, "test_")
testfiles = sort(filter(istest, readdir(testdir)))
if length(testfiles) < 1 error("No tests files found!") end

testfiles = [
    ("test_edge_sync.jl", 4),
    # ("test_indexing.jl", 16),
    # ("test_topology_1d.jl", 16),
    # ("test_topology_2d.jl", 16),
    # ("test_topology_3d.jl", 16),
    # ("test_mpihaloarray.jl", 8),
    ]

# Run each test with the mpiexec command
@testset "$(f[1])" for f in testfiles
    file, nprocs = f
    mpiexec() do cmd
        # e.g. run `mpiexec -n 4 julia test_somthing.jl`
        run(`$cmd -n $nprocs $(Base.julia_cmd()) $(joinpath(testdir, file))`)
        @test true
    end
end

nothing