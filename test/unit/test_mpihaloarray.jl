using MPI, Test
using InteractiveUtils
using MPIHaloArrays


MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)

thisfile = @__FILE__
@assert nprocs == 8 "The tests in $(thisfile) are designed with 8 processes only"

function test_1d_array()
    # Create the MPI topology, which here is a Cartesian 4x4 domain (using 16 cores)
    # topology = CartesianTopology([4,4], [true, true])
    topology = CartesianTopology(8, true)

    # How many halo cells in each dimension (fixed for all dimensions)
    nhalo = 3
    ni = 5

    # Create an uninitialized array with the given topology
    x = MPIHaloArray(1:ni, topology, nhalo)
    if x.topology.rank == 1
        @test all(x.local_indices[1].lo_halo .== (1,3))
        @test all(x.local_indices[1].lo_halo_domain_donor == (4,6))
        @test all(x.local_indices[1].domain == (4,8))
        @test all(x.local_indices[1].hi_halo_domain_donor == (6,8))
        @test all(x.local_indices[1].hi_halo == (9,11))

        @test all(x.global_indices[1].lo_halo .== (12,14))
        @test all(x.global_indices[1].lo_halo_domain_donor == (15,17))
        @test all(x.global_indices[1].domain == (15,19))
        @test all(x.global_indices[1].hi_halo_domain_donor == (17,19))
        @test all(x.global_indices[1].hi_halo == (20,22))
    end
end

function test_2d_array()
    topology = CartesianTopology([4,2], [true,true])

    # How many halo cells in each dimension (fixed for all dimensions)
    nhalo = 2
    ni = 10
    nj = 20

    # Create an uninitialized array with the given topology
    x = MPIHaloArray(rand(Int8,ni,nj), topology, nhalo)

    if x.topology.rank == 1
        # i dim
        @test all(x.topology.coords .== (1,0))

        @test all(x.local_indices[1].lo_halo .== (1,2))
        @test all(x.local_indices[1].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[1].domain == (3,12))
        @test all(x.local_indices[1].hi_halo_domain_donor == (11,12))
        @test all(x.local_indices[1].hi_halo == (13,14))

        @test all(x.global_indices[1].lo_halo .== (15,16))
        @test all(x.global_indices[1].lo_halo_domain_donor == (17,18))
        @test all(x.global_indices[1].domain == (17,26))
        @test all(x.global_indices[1].hi_halo_domain_donor == (25,26))
        @test all(x.global_indices[1].hi_halo == (27,28))

        # j dim
        @test all(x.local_indices[2].lo_halo .== (1,2))
        @test all(x.local_indices[2].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[2].domain == (3,22))
        @test all(x.local_indices[2].hi_halo_domain_donor == (21,22))
        @test all(x.local_indices[2].hi_halo == (23,24))

        # the domain coordinates are (i,j) = (1,0), so no difference
        # for global j indices on this proc
        @test all(x.global_indices[2].lo_halo .== (1,2))
        @test all(x.global_indices[2].lo_halo_domain_donor == (3,4))
        @test all(x.global_indices[2].domain == (3,22))
        @test all(x.global_indices[2].hi_halo_domain_donor == (21,22))
        @test all(x.global_indices[2].hi_halo == (23,24))
    elseif x.topology.rank == 4

        @test all(x.topology.coords .== (0,1))

        # i dim
        @test all(x.local_indices[1].lo_halo .== (1,2))
        @test all(x.local_indices[1].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[1].domain == (3,12))
        @test all(x.local_indices[1].hi_halo_domain_donor == (11,12))
        @test all(x.local_indices[1].hi_halo == (13,14))

        @test all(x.global_indices[1].lo_halo .== (1,2))
        @test all(x.global_indices[1].lo_halo_domain_donor == (3,4))
        @test all(x.global_indices[1].domain == (3,12))
        @test all(x.global_indices[1].hi_halo_domain_donor == (11,12))
        @test all(x.global_indices[1].hi_halo == (13,14))


        # j dim
        @test all(x.local_indices[2].lo_halo .== (1,2))
        @test all(x.local_indices[2].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[2].domain == (3,22))
        @test all(x.local_indices[2].hi_halo_domain_donor == (21,22))
        @test all(x.local_indices[2].hi_halo == (23,24))

        @test all(x.global_indices[2].lo_halo .== (25,26))
        @test all(x.global_indices[2].lo_halo_domain_donor == (27,28))
        @test all(x.global_indices[2].domain == (27,46))
        @test all(x.global_indices[2].hi_halo_domain_donor == (45,46))
        @test all(x.global_indices[2].hi_halo == (47,48))
    end
end

function test_3d_array()
    topology = CartesianTopology([2,2,2], [true,true,true])

    # How many halo cells in each dimension (fixed for all dimensions)
    nhalo = 2
    ni = 10
    nj = 20
    nk = 15

    # Create an uninitialized array with the given topology
    x = MPIHaloArray(rand(Int8,ni,nj,nk), topology, nhalo)

    if x.topology.rank == 1
        # i dim
        @test all(x.topology.coords .== (1,0,0))

        @test all(x.local_indices[1].lo_halo .== (1,2))
        @test all(x.local_indices[1].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[1].domain == (3,12))
        @test all(x.local_indices[1].hi_halo_domain_donor == (11,12))
        @test all(x.local_indices[1].hi_halo == (13,14))

        @test all(x.global_indices[1].lo_halo .== (15,16))
        @test all(x.global_indices[1].lo_halo_domain_donor == (17,18))
        @test all(x.global_indices[1].domain == (17,26))
        @test all(x.global_indices[1].hi_halo_domain_donor == (25,26))
        @test all(x.global_indices[1].hi_halo == (27,28))

        # j dim
        @test all(x.local_indices[2].lo_halo .== (1,2))
        @test all(x.local_indices[2].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[2].domain == (3,22))
        @test all(x.local_indices[2].hi_halo_domain_donor == (21,22))
        @test all(x.local_indices[2].hi_halo == (23,24))

        @test all(x.global_indices[2].lo_halo .== (1,2))
        @test all(x.global_indices[2].lo_halo_domain_donor == (3,4))
        @test all(x.global_indices[2].domain == (3,22))
        @test all(x.global_indices[2].hi_halo_domain_donor == (21,22))
        @test all(x.global_indices[2].hi_halo == (23,24))

        # k dim
        @test all(x.local_indices[3].lo_halo .== (1,2))
        @test all(x.local_indices[3].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[3].domain == (3,17))
        @test all(x.local_indices[3].hi_halo_domain_donor == (16,17))
        @test all(x.local_indices[3].hi_halo == (18,19))

        @test all(x.global_indices[3].lo_halo .== (1,2))
        @test all(x.global_indices[3].lo_halo_domain_donor == (3,4))
        @test all(x.global_indices[3].domain == (3,17))
        @test all(x.global_indices[3].hi_halo_domain_donor == (16,17))
        @test all(x.global_indices[3].hi_halo == (18,19))

    elseif x.topology.rank == 5

        @test all(x.topology.coords .== (1,0,1))

        # i dim
        @test all(x.local_indices[1].lo_halo .== (1,2))
        @test all(x.local_indices[1].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[1].domain == (3,12))
        @test all(x.local_indices[1].hi_halo_domain_donor == (11,12))
        @test all(x.local_indices[1].hi_halo == (13,14))

        @test all(x.global_indices[1].lo_halo .== (15,16))
        @test all(x.global_indices[1].lo_halo_domain_donor == (17,18))
        @test all(x.global_indices[1].domain == (17,26))
        @test all(x.global_indices[1].hi_halo_domain_donor == (25,26))
        @test all(x.global_indices[1].hi_halo == (27,28))

        # j dim
        @test all(x.local_indices[2].lo_halo .== (1,2))
        @test all(x.local_indices[2].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[2].domain == (3,22))
        @test all(x.local_indices[2].hi_halo_domain_donor == (21,22))
        @test all(x.local_indices[2].hi_halo == (23,24))

        @test all(x.global_indices[2].lo_halo .== (1,2))
        @test all(x.global_indices[2].lo_halo_domain_donor == (3,4))
        @test all(x.global_indices[2].domain == (3,22))
        @test all(x.global_indices[2].hi_halo_domain_donor == (21,22))
        @test all(x.global_indices[2].hi_halo == (23,24))

        # k dim
        @test all(x.local_indices[3].lo_halo .== (1,2))
        @test all(x.local_indices[3].lo_halo_domain_donor == (3,4))
        @test all(x.local_indices[3].domain == (3,17))
        @test all(x.local_indices[3].hi_halo_domain_donor == (16,17))
        @test all(x.local_indices[3].hi_halo == (18,19))

        @test all(x.global_indices[3].lo_halo .== (20,21))
        @test all(x.global_indices[3].lo_halo_domain_donor == (22,23))
        @test all(x.global_indices[3].domain == (22,36))
        @test all(x.global_indices[3].hi_halo_domain_donor == (35,36))
        @test all(x.global_indices[3].hi_halo == (37,38))
    end
end

# n halo cells ≥ n real cells
# @test_throws AssertionError MPIHaloArray(zeros(6,6), topology, 6)

# Fill each by the current rank
# fill!(x, rank)

test_1d_array()
test_2d_array()
test_3d_array()

GC.gc()
MPI.Finalize()