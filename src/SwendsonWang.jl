function SWstep!(lattice::T, AD::AdjacencyData, P_nfreeze_fn::Function) where T <: ContinuousLattice
    create_cuts!(lattice)
    update!(AD, lattice)
    freeze_bonds!(lattice, AD, P_nfreeze_fn)
    visit_history = copy(lattice.states)
    
    # Identify clusters and flip
    for (cutidx, spinidx, _) in zip(findnz(lattice.states)...)
        flipfactor = rand() < 0.5 ? -1 : 1
        backtrack(lattice, AD, visit_history, spinidx, cutidx; flipfactor = flipfactor)
    end
    
    clean_cuts!(lattice)
end

function freeze_bonds!(lattice::T, AD::AdjacencyData, P_nfreeze_fn::Function) where {T <: ContinuousLattice}
    for (a, b, idx) in zip(findnz(AD.adjindex_data)...)
        i = min(a, b)
        j = max(a, b)
        mat = AD.adjacency_data[idx]
        
        
        for (i_cut, j_cut, state) in zip(findnz(mat)...)
            if state == -1
                mat[i_cut, j_cut] = 0
                continue
            else
                t = idxtocut(lattice, noverlap(lattice, i, i_cut, j, j_cut))
                if (rand() < P_nfreeze_fn(t))
                    mat[i_cut, j_cut] = 0
                end
            end
        end
    end
    
    for mat in AD.adjacency_data
        dropzeros!(mat)
    end
end


visited(visit_history, spin, cut) = visit_history[cut, spin] == 0

function backtrack(lattice::T, AD::AdjacencyData, visit_history, spin, cut; flipfactor = -1) where T <: ContinuousLattice
    if !visited(visit_history, spin, cut)
        visit_history[cut, spin] = 0
        lattice.states[cut, spin] *= flipfactor
        for (nspin, ncuts) in AD(lattice, spin, cut)
            for ncut in ncuts
                backtrack(lattice, AD, visit_history, nspin, ncut; flipfactor = flipfactor)
            end
        end
    end
end

