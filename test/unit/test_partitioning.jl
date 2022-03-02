include("../../src/partitioning.jl")

@testset "2DHaloPartition_3DArray_2DTiles" begin

    A = zeros(4,100,200);
    tile_dims = (2,2)
    halo_dims = (2,3)
    I, J = get_subdomain_dimension_sizes(size(A), tile_dims, halo_dims)

    @test all(I .== 50)
    @test all(J .== 100)
end

@testset "3DHaloPartition_4DArray_3DTiles" begin

    A = zeros(4,100,200,400);
    tile_dims = (2,2,2)
    halo_dims = (2,3,4)
    I, J, K = get_subdomain_dimension_sizes(size(A), tile_dims, halo_dims)

    @test all(I .== 50)
    @test all(J .== 100)
    @test all(K .== 200)
end

@testset "3DHaloPartition_3DArray_3DTiles" begin

    A = zeros(100,200,400);
    tile_dims = (2,2,2)
    halo_dims = (1,2,3)
    I, J, K = get_subdomain_dimension_sizes(size(A), tile_dims, halo_dims)

    sizes = get_subdomain_sizes(size(A), tile_dims, halo_dims)

    @test all(sizes[1,:,:] .== size(A, 1) / tile_dims[1])
    @test all(sizes[2,:,:] .== size(A, 2) / tile_dims[2])
    @test all(sizes[3,:,:] .== size(A, 3) / tile_dims[3])
    @test all(I .== size(A, 1) / tile_dims[1])
    @test all(J .== size(A, 2) / tile_dims[2])
    @test all(K .== size(A, 3) / tile_dims[3])
end

@testset "3DHaloPartition_3DArray_2DTiles" begin

    A = zeros(100,200,400);
    tile_dims = (2,4,1)
    halo_dims = (1,2,3)
    I, J, K = get_subdomain_dimension_sizes(size(A), tile_dims, halo_dims)
    
    sizes = get_subdomain_sizes(size(A), tile_dims, halo_dims)

    @test all(I .== 50)
    @test all(J .== 50)
    @test all(K .== 400)
end


