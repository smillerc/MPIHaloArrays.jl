
struct DataIndices{T <: Int64}
    lo_halo::NTuple{2,T}               # (start,end) indices of the halo region on the low side
    lo_halo_domain_donor::NTuple{2,T}  # (start,end) indices of the domain region on the low side that "donates" to the halo region
    domain::NTuple{2,T}                # (start,end) indices of the real domain
    hi_halo_domain_donor::NTuple{2,T}  # (start,end) indices of the domain region on the high side that "donates" to the halo region
    hi_halo::NTuple{2,T}               # (start,end) indices of the halo region on the high side
end

localindices(DI::DataIndices) = DI.domain
lohaloindices(DI::DataIndices) = DI.lo_halo
lodonorindices(DI::DataIndices) = DI.lo_halo_domain_donor
hihaloindices(DI::DataIndices) = DI.hi_halo
hidonorindices(DI::DataIndices) = DI.hi_halo_domain_donor
