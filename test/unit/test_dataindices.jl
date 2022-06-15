using Test

include("../../src/MPIHaloArrays.jl")
include("../../src/utils/indexing.jl")
include("../../src/utils/dataindices.jl")

@testset "DataIndices" begin
  N = 2
  A_with_halo = zeros(9, 7)
  nhalo = 2

  local_di = Vector{DataIndices}(undef, N)

  for dim in 1:N
    lo_halo_start, lo_halo_end, lo_dom_start, lo_dom_end = lo_indices(A_with_halo, dim, nhalo)
    hi_dom_start, hi_dom_end, hi_halo_start, hi_halo_end = hi_indices(A_with_halo, dim, nhalo)

    local_di[dim] = DataIndices((lo_halo_start, lo_halo_end),
      (lo_dom_start, lo_dom_end),
      (lo_dom_start, hi_dom_end),
      (hi_dom_start, hi_dom_end),
      (hi_halo_start, hi_halo_end))
  end

  @test all(lohaloindices(local_di[1]) .== (1, 2))
  @test all(lodonorindices(local_di[1]) .== (3, 4))
  @test all(hidonorindices(local_di[1]) .== (6, 7))
  @test all(hihaloindices(local_di[1]) .== (8, 9))
  @test all(domainindices(local_di[1]) .== (3, 7))

  @test all(lohaloindices(local_di[2]) .== (1, 2))
  @test all(lodonorindices(local_di[2]) .== (3, 4))
  @test all(hidonorindices(local_di[2]) .== (4, 5))
  @test all(hihaloindices(local_di[2]) .== (6, 7))
  @test all(domainindices(local_di[2]) .== (3, 5))
end