function setallzero!(x::SparseMatrixCSC{T, I}) where {T, I}
    a, b, _ = findnz(x)
    x[a,b] .= zero(T)
end
function empty!(x::SparseMatrixCSC{T, I}) where {T, I}
    setallzero!(x)
    dropzeros!(x)
end
function lastcutidx(array, cutidx)
    idxs = array
    if cutidx == last(array)
        return cutidx
    end
    for i in eachindex(idxs)
        if idxs[i] > cutidx
            return idxs[i-1]
        end
    end
    return
end

function nextcutidx(array, cutidx)
    for i in eachindex(array)
        if array[i] > cutidx
            return array[i]
        end
    end
    return
end
function isoverlaping(a1, a2, b1, b2)
    if b1 > a2 || a1 > b2
        true
    else
        false
    end
end
function noverlap(a1, a2, b1, b2)
    if b1 > a2 || a1 > b2
        return 0
    else
        return min(a2, b2) - max(a1, b1)
    end
end

