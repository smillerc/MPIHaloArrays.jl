using Test
using MPI

MPI.Init()
comm_size = MPI.Comm_size(MPI.COMM_WORLD)
comm = MPI.COMM_WORLD

function testme(comm)
    print("Hello world, I am rank $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))\n")
end

testme(comm)

MPI.Barrier(MPI.COMM_WORLD)

GC.gc()
MPI.Finalize()

@test MPI.Finalized()
