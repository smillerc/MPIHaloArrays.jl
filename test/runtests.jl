using Test, MPI

nprocs_str = get(ENV, "JULIA_MPI_TEST_NPROCS","")


testdir = @__DIR__
istest(f) = endswith(f, ".jl") && startswith(f, "test_")
testfiles = sort(filter(istest, readdir(testdir)))
if length(testfiles) < 1 error("No tests files found!") end



# Run each test with the mpiexec command
@testset "$f" for f in testfiles
    # if f == "test_topology_2d.jl" || f == "test_topology_3d.jl"
        nprocs = 16
    # else
    #     nprocs = nprocs_str == "" ? clamp(Sys.CPU_THREADS, 2, 4) : parse(Int, nprocs_str) 
    # end
    mpiexec() do cmd
        # e.g. run `mpiexec -n 4 julia test_somthing.jl`
        run(`$cmd -n $nprocs $(Base.julia_cmd()) $(joinpath(testdir, f))`)
        @test true
    end
end

nothing