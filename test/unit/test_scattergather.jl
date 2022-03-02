using MPI, Test
using InteractiveUtils

using MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

# @assert nprocs == 8 "This MPIHaloArray edge_sync test is designed with 6 processes only"

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
    A_local = scatterglobal(A_global, root, nhalo, topology)

    A_global_result = gatherglobal(A_local; root=root)
    if rank == root
        @test all(A_global_result .== A_global)
    end
end

function test_2darray_scatter_gather_detail()
    topology = CartesianTopology(comm, [2,2], [false, false])

    root = 0
    nhalo = 1
    ni = 123 # use random odd numbers for weird domain shapes
    nj = 437 # use random odd numbers for weird domain shapes

    A_global = reshape(1:ni*nj|> collect, ni, nj);
    # if rank == root
    #     display(A_global)
    #     println()
    # end
    A_local = scatterglobal(A_global, root, nhalo, topology)
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
    A_local = scatterglobal(A_global, root, nhalo, topology)

    A_global_result = gatherglobal(A_local; root=root)
    if rank == root
        @test all(A_global_result .== A_global)
    end
end

test_2darray_scatter_gather_detail()
# test_1darray_scatter_gather()
# test_2darray_scatter_gather()
# test_3darray_scatter_gather()

GC.gc()
MPI.Finalize()
