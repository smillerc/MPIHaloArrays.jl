
"""
Perform a global sum operation

# Arguments
 - `A`: MPIHaloArray to perform the operation on
 - `broadcast`: true/false - broadcast to all MPI ranks [default is false]
 - `root`: If `broadcast` is false, which MPI rank to reduce to
"""
function globalsum(A::MPIHaloArray{T,N}; root=0, broadcast=false) where {T,N}
    localsum = sum(domainview(A))

    if broadcast
        globalsum = MPI.Allreduce(localsum, MPI.MPI_SUM, A.topology.comm)
    else
        globalsum = MPI.Reduce(localsum, MPI.MPI_SUM, root, A.topology.comm)
    end

    globalsum
end

"""
Perform a global maximum operation

# Arguments
 - `A`: MPIHaloArray to perform the operation on
 - `broadcast`: true/false - broadcast to all MPI ranks [default is false]
 - `root`: If `broadcast` is false, which MPI rank to reduce to
"""
function globalmax(A::MPIHaloArray{T,N}; root=0, broadcast=false) where {T,N}
    localmax = maximum(domainview(A))

    if broadcast
        globalmax = MPI.Allreduce(localmax, MPI.MPI_MAX, A.topology.comm)
    else
        globalmax = MPI.Reduce(localmax, MPI.MPI_MAX, root, A.topology.comm)
    end
    globalmax
end

"""
Perform a global minimum operation

# Arguments
 - `A`: MPIHaloArray to perform the operation on
 - `broadcast`: true/false - broadcast to all MPI ranks [default is false]
 - `root`: If `broadcast` is false, which MPI rank to reduce to
"""
function globalmin(A::MPIHaloArray{T,N}; root=0, broadcast=false) where {T,N}
    localmin = minimum(domainview(A))

    if broadcast
        globalmin = MPI.Allreduce(localmin, MPI.MPI_MIN, A.topology.comm)
    else
        globalmin = MPI.Reduce(localmin, MPI.MPI_MIN, root, A.topology.comm)
    end
    globalmin
end