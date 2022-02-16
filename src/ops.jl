
function globalsum(A::MPIHaloArray{T,N}; broadcast=false) where {T,N}
    localsum = sum(A.data)
    
    if broadcast
        globalsum = MPI.AllReduce(localsum, +, A.comm)
    else
        globalsum = MPI.Reduce(localsum, +, 0, A.comm)
        MPI.Bcast!(globalsum, root, A.comm)
    end

    globalsum
end

# function Base.maximum(A::AbstractArray; broadcast=false)
#     localmaximum = maximum(A)
# end
# function globalmin(A::MPIHaloArray{T,N}; broadcast=false) where {T,N}
#     localminimum = minimum(A.data)
    
#     globalminimum = MPI.Reduce(localminimum, minimum, 0, A.comm)

#     if broadcast
#     else
#         MPI.Bcast!(globalminimum, root, A.comm)
#     end

#     globalminimum
# end

# function globalmax(A::MPIHaloArray{T,N}; broadcast=false) where {T,N}
#     localmaximum = maximum(A.data)
    
#     globalmaximum = MPI.Reduce(localmaximum, +, 0, A.comm)

#     if broadcast
#         MPI.Bcast!(globalmaximum, root, A.comm)
#     end

#     globalmaximum
# end