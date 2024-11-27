# partitions for multithreading are defined here

"""
    get_partition(M, nthreads)

Return `nthreads` vectors of indices partitioning `M` grid points for multithreading

"""
function get_partition(M, nthreads)

    elements_per_part = fld(M, nthreads)
    partition = Vector{UnitRange{Int64}}(undef, nthreads)
    for id = 1:(nthreads-1)
        partition[id] = (1+elements_per_part*(id-1)):(elements_per_part*id)
    end
    partition[nthreads] = (1+elements_per_part*(nthreads-1)):M
    @assert length(partition) == nthreads
    return partition
end

"""
    get_partition(Mx, My, nthreads, nthx=0, nthy=0)

Return `nthreads` vectors of indices partitioning `Mx` × `My` grid points for multithreading. 

* `nthx`: number of cuts in the ``x`` direction
* `nthy`: number of cuts in the ``y`` direction

We must have `nthx` × `nthy` == `nthreads`.

"""
function get_partition(Mx, My, nthreads, nthx=0, nthy=0)
    if nthx <=0 || nthy<=0
        nthx = gcdx(nthreads,Mx)[1]
        nthy = Int(nthreads/nthx)
    else
        @assert nthx*nthy == nthreads "partition grid size of $nthx x $nthy does not match total number of threads $nthreads"
    end

    partx = get_partition(Mx, nthx)
    party = get_partition(My, nthy)
    partition = Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}}(undef, nthreads)
    for idy = 1:nthy
        for idx = 1:nthx
            id = idx + (idy - 1) * nthx
            partition[id] = (partx[idx], party[idy])
        end
    end
    @assert length(partition) == nthreads
    return partition
end
