mutable struct AdjacencyData
    adjindex_data::SparseMatrixCSC{Int32}
    adjacency_data::Vector{SparseMatrixCSC{Int32}}
end


function AdjacencyData(lattice::T) where {T <: ContinuousLattice}
    nSpins = size(lattice.states, 2)
    maxnseg = round(Int, lattice.beta / lattice.prec)

    adjindex_data  = spzeros(Int32, nSpins, nSpins)
    adjacency_data = SparseMatrixCSC{Int8}[]
    
    for myidx in eachindex(eachcol(lattice.states))
        nzcuts = findnz(lattice.states[:, myidx])[1][1:end-1]

        for cutidx in nzcuts
            mystate = lattice.states[cutidx, myidx]
            for nidx in neighbours(lattice, myidx)
                if myidx > nidx
                    continue
                end
                if iszero(adjindex_data[myidx, nidx])
                    push!(adjacency_data, spzeros(Int32, maxnseg, maxnseg))
                    adjindex_data[myidx, nidx] = lastindex(adjacency_data)
                    mynmat = last(adjacency_data)
                else
                    mynmat = adjacency_data[adjindex_data[myidx, nidx]]
                end
                
                mynzidxs, _ = findnz(lattice.states[:, myidx])
                mycutend_idx = nextcutidx(mynzidxs, cutidx)
                narr = lattice.states[:, nidx]
                narr_nzidxs, _ = findnz(narr)
                ncutidx = lastcutidx(narr_nzidxs, cutidx)
                while !isnothing(ncutidx) && ncutidx < mycutend_idx
                    mynmat[cutidx, ncutidx] = mystate == lattice.states[ncutidx, nidx] ? 1 : -1
                    ncutidx = nextcutidx(narr_nzidxs, ncutidx)
                end
            end
        end
    end

    AdjacencyData(adjindex_data, adjacency_data)
end

function show(::IO, AD::AdjacencyData)
    println("Adjacency data with $(length(AD.adjacency_data)) bonds.")
    return
end

function (AD::AdjacencyData)(chain::T, spinidx, cutidx; aligned = true) where T <: ContinuousLattice
    d = Dict{Int32, Vector{Int32}}()
    if aligned
        eqv = 1
    else
        eqv = -1
    end
    
    for n in neighbours(chain, spinidx)
        if spinidx < n
            mat = AD.adjacency_data[AD.adjindex_data[spinidx, n]]
        else
            mat = AD.adjacency_data[AD.adjindex_data[n, spinidx]]'
        end

        idxs, alignment = findnz(mat[cutidx, :])
        d[n] = idxs[findall(==(eqv), alignment)]
    end
    return d
end


function update!(AD::AdjacencyData, lattice::T) where {T <: ContinuousLattice}
    maxnseg = round(Int, lattice.beta / lattice.prec)


    @assert size(AD.adjacency_data[1]) == (maxnseg, maxnseg) "Incorrect number of dimensions!"

    for myidx in eachindex(eachcol(lattice.states)) # Pick a spin
        N = neighbours(lattice, myidx)
        myarr = lattice.states[:, myidx]
        myarr_nzidxs, _ = findnz(myarr)
        
        for nidx in N
            if myidx > nidx
                continue
            end

            mynmat = AD.adjacency_data[AD.adjindex_data[myidx, nidx]]
            narr = lattice.states[:, nidx]
            narr_nzidxs, _ = findnz(narr)
            setallzero!(mynmat)
            
            for mycutidx in myarr_nzidxs[1:end-1]
                mystate = myarr[mycutidx]
                mycutend_idx = nextcutidx(myarr_nzidxs, mycutidx)
                ncutidx = lastcutidx(narr_nzidxs, mycutidx)
                
                while !isnothing(ncutidx) && ncutidx < mycutend_idx
                    mynmat[mycutidx, ncutidx] = mystate == lattice.states[ncutidx, nidx] ? 1 : -1
                    ncutidx = nextcutidx(narr_nzidxs, ncutidx)
                end
            end
        end
    end

    for x in AD.adjacency_data
        dropzeros!(x)
    end
end
