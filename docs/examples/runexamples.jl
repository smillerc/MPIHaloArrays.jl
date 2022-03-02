using MPI

benchmarkdir = @__DIR__

benchmark_files = [
    ("04-diffusion2d.jl", 4),
    ]

# Run each test with the mpiexec command -- this requires special treatment due to MPI
for f in benchmark_files
    file, np = f
    mpiexec() do cmd
        # e.g. run `mpiexec -n 4 julia test_somthing.jl`
        run(`$cmd -n $np $(Base.julia_cmd()) $(joinpath(benchmarkdir, file))`)
    end
end

nothing