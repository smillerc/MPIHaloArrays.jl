using Documenter, MPIHaloArrays
using MPI

EXAMPLES = [
    "1D Halo Example" => "examples/01-halo1d.md",
    "2D Halo Example" => "examples/02-halo2d.md",
    "3D Halo Example" => "examples/03-halo3d.md",
]

examples_md_dir = joinpath(@__DIR__,"src/examples")
isdir(examples_md_dir) || mkdir(examples_md_dir)

for (example_title, example_md) in EXAMPLES
    example_jl = example_md[1:end-2]*"jl"
    @info "Building $example_md"
    open(joinpath(@__DIR__, "src", example_md), "w") do mdfile
        println(mdfile, "# $example_title")
        println(mdfile)
        println(mdfile, "```julia")
        println(mdfile, "# $example_jl")
        println(mdfile, readchomp(joinpath(@__DIR__,example_jl)))
        println(mdfile, "```")
        println(mdfile)

        println(mdfile, "```")
        println(mdfile, "> mpiexecjl -n 8 julia $example_jl")
        cd(@__DIR__) do
            write(mdfile, mpiexec(cmd -> read(`$cmd -n 8 $(Base.julia_cmd()) --project $example_jl`)))
        end
        println(mdfile, "```")
    end
end


DocMeta.setdocmeta!(MPIHaloArrays, :DocTestSetup, :(using MPIHaloArrays); recursive=true)

makedocs(
    source = "src",
    build = "build",
    clean = true,
    doctest = true,
    repo = "https://github.com/smillerc/MPIHaloArrays.jl",
    highlightsig = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    expandfirst = [],
    modules = [MPIHaloArrays],
    pages = [
        "MPIHaloArrays.jl" => "index.md",
        "Examples" => EXAMPLES,
        "Reference" => [
            # "MPIHaloArray" => "mpihaloarray.md",
            # "Topology" => "topology.md",
            "API" => "api.md"
        ]
    ],
    sitename = "MPIHaloArrays.jl"
)

deploydocs(repo = "github.com/smillerc/MPIHaloArrays.jl.git",
           branch = "gh-pages",
           devbranch = "main",)
