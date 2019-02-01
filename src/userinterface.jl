# ====================================================================
#                           USER INTERFACE
# ====================================================================

# ____________________________________________________________________
# Wrapper for user interface
#
struct Wrap{T}
    _w::T
end
wrappee(v::Wrap{T}) where T = getfield(v, :_w)


# ____________________________________________________________________
# Model
#
Model() = Wrap{Model}(Model(nothing))
Model(args::Vararg{Pair{Symbol, T}, N}) where {T<:AbstractComponent, N} =
    addcomp!(Model(), args...)

propertynames(w::Wrap{Model}) = collect(keys(wrappee(w).comp))

getproperty(w::Wrap{Model}, s::Symbol) = get(wrappee(w).comp, s, nothing)

function getindex(w::Wrap{Model}, id::Int=1)
    model = wrappee(w)
    c = model.instruments
    @assert length(c) >= 1 "No domain in model"
    @assert 1 <= id <= length(c) "Invalid index (allowed range: 1 : " * string(length(c)) * ")"
    return Wrap{Instrument}(c[id])
end

propertynames(w::Wrap{Instrument}) = [:domain; wrappee(w).exprnames; wrappee(w).compnames]

function getproperty(w::Wrap{Instrument}, s::Symbol)
    instr = wrappee(w)
    if s == :domain
        ret = instr.domain
        (ndims(ret) == 1)  &&  (return ret[1])
        return ret
    end
    i = findall(s .== instr.compnames)
    if length(i) == 1
        return instr.compevals[i[1]].result
    end
    i = findall(s .== instr.exprnames)
    if length(i) == 1
        return instr.exprevals[i[1]]
    end
    return nothing
end

getparamvalues(w::Wrap{Model}) = getparamvalues(wrappee(w))
setparamvalues!(w::Wrap{Model}, pval::Vector{Float64}) = setparamvalues!(wrappee(w), pval)

# ____________________________________________________________________
# Components
#
function addcomp!(w::Wrap{Model}, args::Vararg{Pair{Symbol, T}, N}) where {T<:AbstractComponent, N}
    model = wrappee(w)
    for c in args
        model.comp[c[1]] = deepcopy(c[2])
        for (pname, param) in getparams(model.comp[c[1]])
            param._private.model = model
            param._private.cname = c[1]
        end
        model.enabled[c[1]] = true
    end
    return w
end

function setfixed!(w::Wrap{Model}, s::Symbol, flag::Bool=true)
    model = wrappee(w)
    model.enabled[s] = !flag
    return w
end

# ____________________________________________________________________
# Model domains
#
function propertynames(p::Parameter)
    out = collect(fieldnames(Parameter))
    out = out[findall(out .!= :_private)]
end

function setproperty!(p::Parameter, s::Symbol, value)
    bkp = p.expr
    setfield!(p, s, convert(typeof(getfield(p, s)), value))
    if (s == :expr)  &&  (p._private.model != nothing)
        try
            recompile!(p._private.model)
        catch err
            show(err)
            setfield!(p, s, bkp)
            recompile!(p._private.model)
        end
    end
    return value
end


# ____________________________________________________________________
# Model domains
#
function add_dom!(w::Wrap{Model}, dom::AbstractDomain)
    model = wrappee(w)
    push!(model.instruments, Instrument(flatten(dom)))
    return length(model.instruments)
end

function add_dom!(w::Wrap{Model}, args...)
    model = wrappee(w)
    push!(model.instruments, Instrument(Domain(args...)))
    return length(model.instruments)
end

rm_dom!(w::Wrap{Model}, id::Int) = deleteat!(wrappee(w).instruments, id)

dom_count(w::Wrap{Model}) = length(wrappee(w).instruments)


# ____________________________________________________________________
# Expressions
#
addexpr!(w::Wrap{Model}, args...; kw...) = _addexpr!(wrappee(w), args...; kw...)

replaceexpr!(w::Wrap{Model}, label::Symbol, expr::Expr) = replaceexpr!(w::Wrap{Model}, 1, label, expr)

function replaceexpr!(w::Wrap{Model}, id::Int, label::Symbol, expr::Expr)
    model = wrappee(w)
    ii = findall(label .== model.instruments[id].exprnames)
    @assert length(ii) == 1 "No expression labelled $label in domain $id"
    model.instruments[id].exprs[ii[1]] = expr
    _recompile!(model, id)
end

setflag!(w::Wrap{Model}, label::Symbol, flag::Bool) = setflag!(w, 1, label, flag)
function setflag!(w::Wrap{Model}, id::Int, label::Symbol, flag::Bool)
    model = wrappee(w)
    ii = findall(label .== model.instruments[id].exprnames)
    @assert length(ii) == 1 "No expression labelled $label in domain $id"
    model.instruments[id].exprcmp[ii[1]] = flag
    recompile!(w)
end


# ____________________________________________________________________
# Evaluate & fit
#
evaluate!(w::Wrap{Model}) = _evaluate!(wrappee(w))

fit!(w::Wrap{Model}, data::AbstractData; kw...) = _fit!(wrappee(w), [data]; kw...)
fit!(w::Wrap{Model}, data::Vector{T}; kw...) where T<:AbstractMeasures = _fit!(wrappee(w), data; kw...)
fit(w::Wrap{Model}, data::AbstractData; kw...) = _fit(wrappee(w), [data]; kw...)
fit(w::Wrap{Model}, data::Vector{T}; kw...) where T<:AbstractMeasures = _fit(wrappee(w), data; kw...)

function propertynames(w::Wrap{FitResult})
    res = wrappee(w)
    out = collect(fieldnames(typeof(res)))
    out = out[findall((out .!= :bestfit)  .&
                      (out .!= :fitter))]
    out = [out; collect(keys(res.bestfit))]
    return out
end

function getproperty(w::Wrap{FitResult}, s::Symbol)
    res = wrappee(w)
    if s in fieldnames(typeof(res))
        return getfield(res, s)
    end
    return Wrap{FitComp}(res.bestfit[s])
end

propertynames(w::Wrap{FitComp}) = collect(keys(wrappee(w).params))

getproperty(w::Wrap{FitComp}, s::Symbol) = wrappee(w).params[s]


# ____________________________________________________________________
# Miscellaneous
#
function test_component(comp::AbstractComponent, args...; iter=1)
    model = Model(:test => comp)
    add_dom!(model, args...)
    addexpr!(model, :(+test))

    printstyled(color=:magenta, bold=true, "First evaluation:\n")
    @time result = evaluate!(model)
    if iter > 0
        println()
        printstyled(color=:magenta, bold=true, "Further evaluations ($iter):\n")

        @time begin
            for i in 1:iter
                result = evaluate!(model)
                wrappee(model).instruments[1].compevals[1].lastParams[1] = NaN  # Force re-calculation
            end
        end
    end
    show(model)
    return nothing
end
test_component(domain::AbstractCartesianDomain, comp::AbstractComponent, iter=1) =
    test_component(flatten(domain), comp, iter)


code(w::Wrap{Instrument}) = println(wrappee(w).code)


probe(data::AbstractMeasures, args...) = probe([data], args...)

function probe(data::Vector{T}, args::Vararg{Tuple{Parameter, Number},N}; nstep=20) where {T<:AbstractMeasures, N}
    for arg in args
        @assert arg[1]._private.model == args[1][1]._private.model
    end
    model = args[1][1]._private.model
    @assert model != nothing
    @assert length(model.instruments) >= 1
    params = collect(values(getparams(model)))
    pvalues = getfield.(params, :val)

    if isa(nstep, Int)
        cd = CartesianDomain(fill(nstep, N)...)
    else
        cd = CartesianDomain(nstep...)
    end
    @assert length(size(cd)) == N

    ipar = Vector{Int}()
    step = Vector{Float64}()
    for arg in args
        for ii in 1:length(params)
            if params[ii] == arg[1]
                push!(ipar, ii)
                push!(step, 2 * arg[2])
                break
            end
        end
    end
    @assert length(ipar) == N
    
    for ii in 1:N
        cd[ii] .-= length(cd[ii])/2. .+ 0.5
        cd[ii] .*= step[ii] / (length(cd[ii])-1)
        cd[ii] .+= params[ipar[ii]].val
    end
    dom = flatten(cd)

    # Inner function to evaluate all the models and store the result in a 1D array
    data1d = data1D(model, data)
    function eval_residuals1d(pvalues::Vector{Float64},
                              model=model, data1d=data1d)
        _evaluate!(model, pvalues)
        return residuals1d(model, data1d)
    end

    out = Matrix{Float64}(undef, length(dom), N+1)
    for ii in 1:length(dom)
        pvalues[ipar] .= dom.axis[:,ii]
        eval_residuals1d(pvalues)
        out[ii, 1:N] .= dom.axis[:,ii]
        out[ii, N+1] = sum(abs2, model.buffer1d)
    end
    return out
end
