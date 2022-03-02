using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

@assert nprocs == 8 "This MPIHaloArray edge_sync test is designed with 6 processes only"

function print_array(U,proc)
    MPI.Barrier(comm)
    if rank == proc
        println("Proc: $proc")
        # for j in 1:size(U.data,2)
        #     println(U.data[:,j])
        # end
        display(U.data)
        println()
    end
    MPI.Barrier(comm)
end

function test_1darray_scatter_gather()
    topology = CartesianTopology(comm, 8, false)

    root = 0
    nhalo = 3
    ni = 259 # use random odd numbers for weird domain shapes

    A_global = 1:ni;
    A_local = scatterglobal(A_global, root, nhalo, topology)

    A_global_result = gatherglobal(A_local; root=root)
    if rank == root
        @test all(A_global_result .== A_global)
    end
end

function test_2darray_scatter_gather()
    topology = CartesianTopology(comm, [4,2], [false, false])

    root = 0
    nhalo = 3
    ni = 293 # use random odd numbers for weird domain shapes
    nj = 131 # use random odd numbers for weird domain shapes

    A_global = reshape(1:ni*nj, ni, nj);

    rank_0_view = @view A_global[1:74, 1:66]
    rank_6_view = @view A_global[148:220, 67:131]

    A_local = scatterglobal(A_global, root, nhalo, topology)

    ilo,ihi,jlo,jhi = localindices(A_local)

    # Test to see if the splitting/partitioning was done correctly
    if rank == 0
        @test all(A_local[ilo:ihi,jlo:jhi] .== rank_0_view)
    elseif rank == 6
        @test all(A_local[ilo:ihi,jlo:jhi] .== rank_6_view)
    end

    A_global_result = gatherglobal(A_local; root=root)
    if rank == root
        @test all(A_global_result .== A_global)
    end
end

function test_3darray_scatter_gather()
    topology = CartesianTopology(comm, [2,2,2], [false, false, false])

    root = 0
    nhalo = 3
    ni = 293 # use random odd numbers for weird domain shapes
    nj = 131 # use random odd numbers for weird domain shapes
    nk = 36  # use random odd numbers for weird domain shapes

    A_global = reshape(1:ni*nj*nk, (ni, nj, nk));

    rank_0_view = @view A_global[1:147, 1:66, 1:18]
    rank_7_view = @view A_global[148:293, 67:131, 19:36]

    A_local = scatterglobal(A_global, root, nhalo, topology)

    ilo, ihi, jlo, jhi, klo, khi = localindices(A_local)

    # Test to see if the splitting/partitioning was done correctly
    if rank == 0
        @test all(A_local[ilo:ihi,jlo:jhi,klo:khi] .== rank_0_view)
    elseif rank == 7
        @test all(A_local[ilo:ihi,jlo:jhi,klo:khi] .== rank_7_view)
    end

    A_global_result = gatherglobal(A_local; root=root)
    if rank == root
        @test all(A_global_result .== A_global)
    end
end

test_1darray_scatter_gather()
test_2darray_scatter_gather()
test_3darray_scatter_gather()

GC.gc()
MPI.Finalize()
