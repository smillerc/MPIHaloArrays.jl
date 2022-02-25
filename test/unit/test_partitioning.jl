include("../../src/partitioning.jl")

@testset "Partition1D" begin
    dims = 13
    ntiles = 4
    indices = [
        (1,3),
        (4,6),
        (7,9),
        (10,13),
    ]
    for id in 1:ntiles
        ilo, ihi = tile_indices_1d(dims, ntiles, id)
        @test all(indices[id] .== (ilo, ihi))
    end
end

@testset "Partition2D" begin
    dims = (28,38)
    ntiles = 6
    indices = [
       (1, 9, 1, 19),
       (10, 18, 1, 19),
       (19, 28, 1, 19),
       (1, 9, 20, 38),
       (10, 18, 20, 38),
       (19, 28, 20, 38),
    ]
    for id in 1:ntiles
        ilo, ihi, jlo, jhi = tile_indices_2d(dims, ntiles, id)
        @test all(indices[id] .== (ilo, ihi, jlo, jhi))
    end
end


@testset "Partition3D" begin
    dims = (28,38,16)
    ntiles = 6
    indices = [
       (1, 9, 1, 19, 1, 16),
       (10, 18, 1, 19, 1, 16),
       (19, 28, 1, 19, 1, 16),
       (1, 9, 20, 38, 1, 16),
       (10, 18, 20, 38, 1, 16),
       (19, 28, 20, 38, 1, 16),
    ]
    for id in 1:ntiles
        ilo, ihi, jlo, jhi, klo, khi = tile_indices_3d(dims, ntiles, id)
        @test all(indices[id] .== (ilo, ihi, jlo, jhi, klo, khi))
    end
end