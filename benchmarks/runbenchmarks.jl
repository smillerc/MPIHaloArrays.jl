using Test, MPI, BenchmarkTools

benchmarkdir = @__DIR__

benchmark_files = [
    ("2d_halo_exchange.jl", 8),
    ]

# Run each test with the mpiexec command -- this requires special treatment due to MPI
@testset "$(f[1])" for f in benchmark_files
    file, nprocs = f
    mpiexec() do cmd
        # e.g. run `mpiexec -n 4 julia test_somthing.jl`
        run(`$cmd -n $nprocs $(Base.julia_cmd()) $(joinpath(benchmarkdir, file))`)
        @test true
    end
end

nothing