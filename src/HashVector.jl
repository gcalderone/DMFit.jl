struct HashVector{T}  <: AbstractDict{Symbol, T}
    keys::Vector{Symbol}
    values::Vector{T}
end

length(v::HashVector{T}) where T = length(v.values)

function HashVector{T}() where T
    return HashVector{T}(Vector{Symbol}(), Vector{T}())
end

function keyToIndex(v::HashVector{T}, key::Symbol) where T
    i = findall(v.keys .== key)
    @assert length(i) <= 1
    if length(i) == 1
        return i[1]
    end
    return 0
end

function push!(v::HashVector{T}, key::Symbol, value::T) where T
    @assert keyToIndex(v, key) == 0 "Key :$key already exists"
    push!(v.keys, key)
    push!(v.values, value)
    return v
end

function HashVector{T}(p::Pair{Symbol,T}) where T
    ret = HashVector{T}(Vector{Symbol}(), Vector{T}())
    push!(ret, p[1], p[2])
    return ret
end

function HashVector{T}(t::Pair{Symbol,T}...) where T
    ret = HashVector{T}(Vector{Symbol}(), Vector{T}())
    for p in t
        push!(ret, p[1], p[2])
    end
    return ret
end

keys(v::HashVector{T}) where T = return v.keys
values(v::HashVector{T}) where T = return v.values

setindex!(v::HashVector{T}, value::T, i::Int) where T = v.values[i] = value
function setindex!(v::HashVector{T}, value::T, key::Symbol) where T
    i = keyToIndex(v, key)
    @assert i != 0 "Can't add data to a HashVector structure via [] syntax"
    v.values[i] = value
end

getindex(v::HashVector{T}, i::Int) where T = v.values[i]
function getindex(v::HashVector{T}, key::Symbol) where T
    i = keyToIndex(v, key)
    @assert i != 0 "Key :$key is not present"
    return v.values[i]
end


if oldver()
    start(v::HashVector{T}) where T = 1
    next(v::HashVector{T}, i::Int) where T = (Pair(v.keys[i], v.values[i]), i + 1)
    done(v::HashVector{T}, i::Int) where T = (i > length(v.values))
else
    function iterate(v::HashVector{T}) where T
        (length(v.values) == 0)  &&  (return nothing)
        i = 1
        return (Pair(v.keys[i], v.values[i]), i)
    end

    function iterate(v::HashVector{T}, i::Int) where T
        (length(v.values) <= i)  &&  (return nothing)
        i += 1
        return (Pair(v.keys[i], v.values[i]), i)
    end
end

#=
if false
    a = ModelFit.HashVector{Float64}()
    push!(a, :a, 1.2)
    push!(a, :b, 1.2)
    a[:a]
    a[:b]
    a[:a] = 34.5
    a[:a]
    a[1] = 12.
    a[:a]
    a[1]

    a = ModelFit.HashVector{Float64}(:a => 1.3)
    a[:a]

    a = ModelFit.HashVector{Float64}(:a => 1.3, :b=>2.4, :c=>3.6)
    a[:a]
    a

    for (sym, val) in a
        println(sym, " ", val)
    end
end
=#
