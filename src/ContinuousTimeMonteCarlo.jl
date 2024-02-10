module ContinuousTimeMonteCarlo

using StatsBase
using SparseArrays
using Memoize
using LazyArrays
using LRUCache
import Base: show, empty!

export AdjacencyData, update!, SWstep!, GenericIsingLattice

abstract type ContinuousLattice end

include("Utilities.jl")
include("Adjacency.jl")
include("SwendsonWang.jl")

function energy(lattice::T, AD::AdjacencyData; refresh = true) where T <: ContinuousLattice
    if refresh
        update!(AD, lattice)
    end
    
    energy = 0
    for (i, j, matidx) in zip(findnz(AD.adjindex_data)...)
        mat = AD.adjacency_data[matidx]
        for (i_cut, j_cut, bondstate) in zip(findnz(mat)...)
            energy -= bondstate * idxtocut(lattice, noverlap(lattice, i, i_cut, j, j_cut))
        end
    end
    energy / maximum(lattice.trange)
end

function magnetization(lattice::T, spinidx) where T <: ContinuousLattice
    segments = lattice.states[:, spinidx]
    total = 0
    for (cutidx, state) in zip(findnz(segments)...)
        seglen = cut_to_range(lattice, spinidx, cutidx)
        total += state * seglen
    end
    total
end
function magnetization(lattice::T) where T <: ContinuousLattice
    total = 0
    for spinidx in eachindex(eachcol(lattice.states))
        total += magnetization(lattice, spinidx)
    end
    total / maximum(lattice.trange)
end

function cut_to_range(lattice::T, spinidx, cutidx) where T <: ContinuousLattice
    states = lattice.states[:, spinidx]
    cutendidx = nextcutidx(findnz(states)[1], cutidx)
    if isnothing(cutendidx)
        return 0
    end
    idxtocut(lattice, (cutendidx - cutidx))
end

function addcut!(chain::T, spinidx, cut) where T <: ContinuousLattice
    states = chain.states[:, spinidx]
    idx = cuttoidx(chain, cut)

    if states[idx] != 0
        return
    end


    nzlst = findnz(states)[1]
    lastidx = lastcutidx(nzlst, idx)
    if isnothing(lastidx)
        @show spinidx, idx
    end

    state = states[lastidx]

    chain.states[idx, spinidx] = Int8(state)
end

idxtocut(lattice::T, index) where {T <: ContinuousLattice} = round(index * lattice.prec, digits = lattice.digits)
cuttoidx(lattice::T, index) where {T <: ContinuousLattice} = round(Int, index / lattice.prec)

function noverlap(chain::T, s1, a1, s2, b1) where {T <: ContinuousLattice}
    nza = findnz(chain.states[:, s1])[1]
    nzb = findnz(chain.states[:, s2])[1]
    a2i = findfirst(==(a1), nza)
    b2i = findfirst(==(b1), nzb)

    if isnothing(a2i) || isnothing(b2i)
        @show isnothing(a2i)
        @show isnothing(b2i)
        @show s1, a1, s2, b1
        return nza, nzb
        throw(Exception("Couldn't find cut index"))
    end
    a2 = nza[a2i+1]
    b2 = nzb[b2i+1]

    noverlap(a1, a2, b1, b2)
end

function create_cuts!(chain::T) where T <: ContinuousLattice
    for (spinidx, _) in enumerate(eachcol(chain.states))
        cut_start = chain.prec
        cut_end = cut_start + sample(chain.trange, chain.w)
        while cut_end < chain.beta
            addcut!(chain, spinidx, cut_end)
            cut_start = cut_end
            cut_end = cut_start + sample(chain.trange, chain.w)
        end
    end
end

function clean_cuts!(chain::T) where {T <: ContinuousLattice}
    for states in eachcol(chain.states)
        nzarr = findnz(states)[1]
    
        nzremoved = 0
        for i in 2:length(nzarr)-1
            curidx = nzarr[i]
            preidx = nzarr[i-(1+nzremoved)]
            
            if states[preidx] == states[curidx]
                states[curidx] = Int8(0)
                nzremoved += 1
            else
                nzremoved = 0
            end
        end
    end    
    dropzeros!(chain.states)
end

function neighbours(lattice::T, linidx) where {T <: ContinuousLattice}
    append!(findnz(lattice.adjmat[linidx, :])[1], findnz(lattice.adjmat[:, linidx])[1])
end

mutable struct GenericIsingLattice{trangeType, weightType} <: ContinuousLattice
    states::SparseMatrixCSC{Int8}
    beta::Float64
    h::Float64
    trange::trangeType
    w::weightType
    prec::Float64
    digits::Int64
    adjmat::SparseMatrixCSC{Int8}
    
    function GenericIsingLattice(N, beta, magnetic_field, adjmat; prec = 1e-3)
        @assert size(adjmat) ==  (N, N) "Adjacency Matrix of incorrect size!"
        states = spzeros(Int8, ceil(Int, beta / prec), N)
        @memoize LRU{Tuple{Float64}, Float64}(maxsize=32) P_timeadd(t) = exp(-magnetic_field * t)
        t = 0:prec:beta
        w = Weights(BroadcastVector(P_timeadd, t))
        
        choices = [Int8(-1), Int8(1)]
        for (_, spinstate) in enumerate(eachcol(states))
            spinstate[1] = rand(choices)
            cut_start = prec
            cut_end = cut_start + sample(t, w)
            while cut_end < beta
                cutidx = round(Int, cut_end / prec)
                spinstate[cutidx] = rand(choices)
                cut_start = cut_end
                cut_end = cut_start + sample(t, w)
            end
            spinstate[end] = Int8(1)
        end
         
        new{typeof(t), typeof(w)}(states, beta, magnetic_field, t, w, prec, Int(log10(prec) |> abs), adjmat)
    end
end

end
