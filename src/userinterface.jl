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



function probe(lparams::Vector{Parameter}, data::Vector{T}; delta=3., kw...) where T <: AbstractMeasures
    (length(lparams) == 0)  &&  (return nothing)
    (length(data) == 0)  &&  (return nothing)
    for par in lparams
        @assert par._private.model == lparams[1]._private.model
    end
    model = lparams[1]._private.model
    @assert model != nothing
    @assert length(model.instruments) >= 1
    params = collect(values(getparams(model)))
    pvalues0 = getfield.(params, :val)
    pvalues = deepcopy(pvalues0)
    
    data1d = data1D(model, data)
    ipar = Vector{Int}()
    range = Vector{Float64}()
    _evaluate!(model, pvalues0)
    cost0 = sum(abs2, residuals1d(model, data1d))
    for par in lparams
        for ii in 1:length(params)
            if params[ii] == par
                push!(ipar, ii)
                @assert isfinite(par._private.fitunc) "The model has not been fit: can't guess step for parameter"
                pvalues = deepcopy(pvalues0)
                lo = pvalues0[ii] - par._private.fitunc/2
                while true
                    if lo < par.low
                        lo = par.low
                        break
                    end
                    pvalues[ii] = lo
                    _evaluate!(model, pvalues)
                    cost = sum(abs2, residuals1d(model, data1d))
                    (abs(cost - cost0) >= delta)  &&  break
                    lo -= par._private.fitunc/2
                end
                hi = pvalues0[ii] + par._private.fitunc
                while true
                    if lo > par.high
                        hi = par.high
                        break
                    end
                    pvalues[ii] = hi
                    _evaluate!(model, pvalues)
                    cost = sum(abs2, residuals1d(model, data1d))
                    (abs(cost - cost0) >= delta)  &&  break
                    hi += par._private.fitunc/2
                end
                append!(range, [lo, hi])
            end
        end
    end
    @assert length(ipar) == length(lparams)
    _evaluate!(model, pvalues0)
    rr = collect(reshape(range, 2, div(length(range),2))')
    return probe(lparams, data, rr; kw...)
end



function probe(lparams::Vector{Parameter}, data::Vector{T}, rr::Matrix{Float64}; nstep=11) where {T<:AbstractMeasures}
    (length(lparams) == 0)  &&  (return nothing)
    (length(data) == 0)  &&  (return nothing)
    for par in lparams
        @assert par._private.model == lparams[1]._private.model
    end
    model = lparams[1]._private.model
    @assert model != nothing
    @assert length(model.instruments) >= 1
    params = collect(values(getparams(model)))
    pvalues0 = getfield.(params, :val)
    pvalues = deepcopy(pvalues0)   
    @assert size(rr)[1] == length(lparams)

    data1d = data1D(model, data)
    _evaluate!(model, pvalues0)
    cost0 = sum(abs2, residuals1d(model, data1d))
   
    tmp = Vector{AbstractArray}()
    ipar = Vector{Int}()
    for par in lparams
        for ii in 1:length(params)
            if params[ii] == par
                push!(ipar, ii)
                push!(tmp, range(rr[length(ipar), 1], stop=rr[length(ipar), 2], length=nstep))
            end
        end
    end
    @assert length(ipar) == length(lparams)

    if length(lparams) > 1
        cd = CartesianDomain(tmp...)
    else
        cd = Domain(tmp...)
    end
    dom = flatten(cd)
    out = Matrix{Float64}(undef, length(dom), length(lparams)+1)
    for ii in 1:length(dom)
        if length(lparams) > 1
            pvalues[ipar] .= dom.axis[:,ii]
            out[ii, 1:end-1] .= dom.axis[:,ii]
        else
            pvalues[ipar] .= dom.axis[ii]
            out[ii, 1] = dom.axis[ii]
        end
        _evaluate!(model, pvalues)
        out[ii, end] = sum(abs2, residuals1d(model, data1d)) .- cost0
    end
    _evaluate!(model, pvalues0)    
    return out
end
