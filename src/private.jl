# ====================================================================
#                         PRIVATE FUNCTIONS
# ====================================================================

# ____________________________________________________________________
# Model
#
function getparams(model::Model)
    out = OrderedDict{Symbol, WParameter}()
    for (cname, comp) in model.comp
        for (pname, par) in getparams(comp, model.enabled[cname])
            out[Symbol(cname, compsep, pname)] = par
        end
    end
    return out
end

getparamvalues(v::Union{Model, AbstractComponent}) = [wpar.par.val for (pname, wpar) in getparams(v)]
function setparamvalues!(v::Union{Model, AbstractComponent}, pval::Vector{Float64})
    i = 0
    for (pname, wpar) in getparams(v)
        i += 1
        wpar.par.val = pval[i]
    end
end


function prepareindex1D!(model::Model)
    out = Vector{Int}()
    push!(out, 0)
    for id in 1:length(model.instruments)
        for jj in 1:length(model.instruments[id].exprnames)
            (model.instruments[id].exprcmp[jj])  ||  (continue)
            push!(out, out[end] + length(model.instruments[id].exprevals[jj]))
        end
    end
    model.index1d = out
end


function model1D(model::Model, buffer=Vector{Float64}())
    (length(buffer) == 0)  &&  (buffer = Vector{Float64}(undef, model.index1d[end]))
    ii = 1
    for id in 1:length(model.instruments)
        for jj in 1:length(model.instruments[id].exprnames)
            (model.instruments[id].exprcmp[jj])  ||  (continue)
            buffer[model.index1d[ii]+1:model.index1d[ii+1]] .= model.instruments[id].exprevals[jj]
            ii += 1
        end
    end
    return buffer
end


function data1D(model::Model, data::Vector{T}) where T<:AbstractMeasures
    @assert length(data) >= (length(model.index1d)-1) "Not enough dataset(s) for model"
    @assert length(data) <= (length(model.index1d)-1) "Too many dataset(s) for model"

    out = Vector{Measures_1D}()
    ii = 1
    for id in 1:length(model.instruments)
        for jj in 1:length(model.instruments[id].exprnames)
            (model.instruments[id].exprcmp[jj])  ||  (continue)
            tmp = length(model.instruments[id].exprevals[jj])
            @assert(length(data[ii]) == tmp,
                    "Length of dataset $ii do not match corresponding model: " *
                    string(length(data[ii])) * " != " * string(tmp))
            d1 = flatten(data[ii], model.instruments[id].domain)
            if ii == 1
                push!(out, d1)
            else
                append!(out[1], d1)
            end
            ii += 1
        end
    end
    return out[1]
end


# ____________________________________________________________________
# Components
#
function getparams(comp::AbstractComponent, enabled::Bool=true)
    out = OrderedDict{Symbol, WParameter}()
    for pname in fieldnames(typeof(comp))
        par = getfield(comp, pname)
        if typeof(par) == Parameter
            (enabled)  ||  (par.fixed = true)
            if par.log
                par.val  = log10(par.val)
                par.low  = log10(par.low)
                par.high = log10(par.high)
            end
            out[pname] = WParameter(pname, 0, par)
        elseif typeof(par) == Vector{Parameter}
            for i in 1:length(par)
                (enabled)  ||  (par[i].fixed = true)
                out[Symbol(pname, i)] = WParameter(pname, i, par[i])
            end
        end
    end
    return out
end


# ____________________________________________________________________
# Instruments
#
# Instrument constructors
Instrument(dom::AbstractDomain) = 
    Instrument(flatten(dom), dom, "", ()->nothing, 0,
               Vector{Symbol}(), Vector{CompEvaluation}(),
               Vector{Symbol}(), Vector{Expr}(), Vector{Bool}(), Vector{Vector{Float64}}())

function Instrument(domain::AbstractDomain, model::Model,
                    labels::Vector{Symbol}, exprs::Vector{Expr}, cmp::Vector{Bool})
    (length(cmp) == 0)  && (cmp = fill(true, length(exprs)))
    @assert length(labels) == length(exprs)
    @assert length(cmp) == length(exprs)

    function parse_model_expr(expr::Union{Symbol, Expr}, cnames, accum=Vector{Symbol}())
        if typeof(expr) == Expr
            # Parse the expression to check which components are involved
            for i in 1:length(expr.args)
                arg = expr.args[i]

                if typeof(arg) == Symbol
                    if arg in cnames # if it is one of the model components...
                        push!(accum, arg)
                    end
                elseif typeof(arg) == Expr
                    parse_model_expr(arg, cnames, accum)
                end
            end
        else
            push!(accum, expr)
        end
        return accum
    end

    # Check which components are involved
    compinvolved = Vector{Symbol}()
    for expr in exprs
        compinvolved = unique([compinvolved; parse_model_expr(expr, keys(model.comp))])
    end

    # Sort involved components according to insertion order
    kk = keys(model.comp)
    sort!(compinvolved, lt=(a, b) -> findall(a .== kk)[1] < findall(b .== kk)[1])

    # Prepare the code for model evaluation
    code = Vector{String}()
    tmp = ""
    for (cname, comp) in model.comp
        for (pname, wpar) in getparams(comp)
            tmp *= ", $(cname)$(compsep)$(pname)::Float64"
        end
    end
    push!(code, "(_instr::Instrument $tmp, _unused_...) -> begin")
    push!(code, "  domain = _instr.ldomain")
    # The last argument, _unused_, allows to push! further
    # components in the model after an expression has already been
    # prepared
    for (cname, comp) in model.comp
        if cname in compinvolved
            i = findall(cname .== compinvolved)[1]
            tmp = ""
            for (pname, wpar) in getparams(comp)
                par = wpar.par
                (par.expr != "")  &&  (push!(code, "  $(cname)$(compsep)$(pname) = " *
                                             replace(par.expr, "this$(compsep)" => "$(cname)$(compsep)")))
                tmp *= ", $(cname)$(compsep)$(pname)"
            end
            if isfunction(comp)
                push!(code, "  $cname(args...) = _evaluate!(_instr.compevals[$i], domain $tmp, args...)")
            else
                push!(code, "  $cname = _evaluate!(_instr.compevals[$i], domain $tmp)")
            end
        end
    end

    exprevals = Vector{Vector{Float64}}()
    if length(exprs) > 1
        for i in 1:length(exprs)
            push!(code, "  " * string(labels[i]) * " = " * string(exprs[i]))
        end
        for i in 1:length(exprs)
            push!(code, "  if length(_instr.exprevals[$i]) == 0")
            push!(code, "    _instr.exprevals[$i]  = " * string(labels[i]))
            push!(code, "  else")
            push!(code, "    _instr.exprevals[$i] .= " * string(labels[i]))
            push!(code, "  end")
            push!(exprevals, Vector{Float64}(undef, 0))
        end
    else
        i = 1
        push!(code, "  if length(_instr.exprevals[$i]) == 0")
        push!(code, "    _instr.exprevals[$i]  = " * string(exprs[i]))
        push!(code, "  else")
        push!(code, "    _instr.exprevals[$i] .= " * string(exprs[i]))
        push!(code, "  end")
        push!(exprevals, Vector{Float64}(undef, 0))
    end
    push!(code, "  _instr.counter += 1")
    push!(code, "  return nothing")
    push!(code, "end")
    funct = eval(Meta.parse(join(code, "\n")))

    compevals = Vector{CompEvaluation}()
    for (cname, comp) in model.comp
        if cname in compinvolved
            cd = cdata(comp, domain)
            npar = length(getparams(comp))
            tmp = CompEvaluation(true, npar, fill(false, npar), 0, cd,
                                 Vector{Float64}(undef, npar),
                                 Vector{Float64}(undef, outsize(cd, domain)))
            tmp.lastParams .= NaN
            push!(compevals, tmp)
        end
    end

    instr = Instrument(domain)
    instr.code = join(code, "\n")
    instr.funct = funct
    instr.counter = 0
    instr.compnames = compinvolved
    instr.compevals = compevals
    instr.exprnames = deepcopy(labels)
    instr.exprs = deepcopy(exprs)
    instr.exprcmp = cmp
    instr.exprevals = exprevals
    return instr
end

# Recompile an instrument
Instrument(model::Model, instr::Instrument) =
    Instrument(instr.domain, model, instr.exprnames, instr.exprs, instr.exprcmp)

function _recompile!(model::Model, id::Int)
    tmp = Instrument(model, model.instruments[id])
    deleteat!(model.instruments, id)
    insert!(model.instruments, id, tmp)
    _evaluate!(model)
    prepareindex1D!(model)
end

# ____________________________________________________________________
# Expressions
#
newexprlabel(model::Model, id::Int) = Symbol(:expr, length(model.instruments[id].exprs)+1)
newexprlabel(model::Model, id::Int, n::Int) = Symbol.(Ref(:expr), length(model.instruments[id].exprs) .+ collect(1:n))

_addexpr!(model::Model                                 , expr::Expr         ; cmp=true) = _addexpr!(model,  1, [newexprlabel(model, 1)]              , [expr]       ; cmp=[cmp])
_addexpr!(model::Model                                 , symbol::Symbol     ; cmp=true) = _addexpr!(model,  1, [newexprlabel(model, 1)]              , [:(+$symbol)]; cmp=[cmp])
_addexpr!(model::Model                                 , exprs::Vector{Expr}; cmp=true) = _addexpr!(model,  1,  newexprlabel(model, 1, length(exprs)), exprs        ; cmp=fill(cmp, length(exprs)))

_addexpr!(model::Model, id::Int                        , expr::Expr         ; cmp=true) = _addexpr!(model, id, [newexprlabel(model, 1)]              , [expr]       ; cmp=[cmp])
_addexpr!(model::Model, id::Int                        , symbol::Symbol     ; cmp=true) = _addexpr!(model, id, [newexprlabel(model, 1)]              , [:(+$symbol)]; cmp=[cmp])
_addexpr!(model::Model, id::Int                        , exprs::Vector{Expr}; cmp=true) = _addexpr!(model, id,  newexprlabel(model, 1, length(exprs)), exprs        ; cmp=fill(cmp, length(exprs)))

_addexpr!(model::Model,          label::Symbol         , expr::Expr         ; cmp=true) = _addexpr!(model,  1, [label]                               , [expr]       ; cmp=[cmp])
_addexpr!(model::Model,          label::Symbol         , symbol::Symbol     ; cmp=true) = _addexpr!(model,  1, [label]                               , [:(+$symbol)]; cmp=[cmp])
_addexpr!(model::Model,          labels::Vector{Symbol}, exprs::Vector{Expr}; cmp=true) = _addexpr!(model,  1,  labels                               , exprs        ; cmp=fill(cmp, length(exprs)))

function _addexpr!(model::Model, id::Int, labels::Vector{Symbol}, exprs::Vector{Expr}; cmp=Vector{Bool}())
    instr = model.instruments[id]
    append!(instr.exprnames, labels)
    append!(instr.exprs, exprs)
    append!(instr.exprcmp, cmp)
    _recompile!(model, id)
    return labels[end]
end


# ____________________________________________________________________
function _evaluate!(c::CompEvaluation, d::AbstractDomain, args...)
    c.enabled  ||  (return c.result)
    if c.counter == 0
        c.counter += 1
        evaluate!(c.cdata, c.result, d, args...)
    else
        #evaluate!(c.cdata, c.result, d, args...)
        for i in 1:c.npar
            (c.log[i])  &&  (args[i] = 10. ^args[i])
            if c.lastParams[i] != args[i]
                c.lastParams .= args
                c.counter += 1
                evaluate!(c.cdata, c.result, d, args...)
            end
            (c.log[i])  &&  (args[i] = log10(args[i]))
        end
    end
    return c.result
end

_evaluate!(instr::Instrument, pvalues::Vector{Float64}) =
    Base.invokelatest(instr.funct, instr, pvalues...)

function _evaluate!(model::Model, pvalues::Vector{Float64})
    for ii in 1:length(model.instruments)
        _evaluate!(model.instruments[ii], pvalues)
    end
    return model
end

function _evaluate!(model::Model)
    for instr in model.instruments
        for ii in 1:length(instr.compnames)
            instr.compevals[ii].enabled = model.enabled[instr.compnames[ii]]
            comp = model.comp[instr.compnames[ii]]
            jj = 1
            for (pname, par) in getparams(comp)
                instr.compevals[ii].log[jj] = par.par.log
            end
        end
    end
    _evaluate!(model, getparamvalues(model))
end


# ____________________________________________________________________
# Minimizer
#
support_param_limits(f::AbstractMinimizer) = false

function _fit(model::Model, data::Vector{T}; kw...) where T<:AbstractMeasures
    pval = getparamvalues(model)
    res = _fit!(model, data; kw...)
    setparamvalues!(model, pval)
    _evaluate!(model)
    return res
end

function _fit!(model::Model, data::Vector{T}; dry=false, minimizer=Minimizer()) where T<:AbstractMeasures
    elapsedTime = Base.time_ns()

    @assert typeof(minimizer) <: AbstractMinimizer
    @assert length(model.instruments) >= 1

    params = getfield.(values(getparams(model)), :par)
    pvalues = getfield.(params, :val)
    uncert = fill(NaN, length(pvalues))
    ifree = findall(.! getfield.(params, :fixed))
    @assert length(ifree) > 0 "No free parameter in the model"

    # Prepare 1D arrays containing all the data and model results
    data1d = data1D(model, data)
    model1d = model1D(model)
    Rmodel1d = Ref(model1d)
    
    # Inner function to evaluate all the models and store the result in a 1D array
    function evaluate1D(freepvalues::Vector{Float64})
        pvalues[ifree] .= freepvalues
        _evaluate!(model, pvalues)
        model1D(model, Rmodel1d[])
        return Rmodel1d[]
    end

    status = :NonOptimal
    if !dry
        # Check if the minimizer supports bounded parameters
        if !support_param_limits(minimizer)
            if  (length(findall(isfinite.(getfield.(params, :low )))) > 0)  ||
                (length(findall(isfinite.(getfield.(params, :high)))) > 0)
                @warn "Parameter bounds are not supported by " * string(typeof(minimizer))
            end
        end

        #Main.code_warntype(evaluate1D, (Vector{Float64},))
        (status, bestfit_val, bestfit_unc) = minimize(minimizer, evaluate1D, data1d.val, data1d.unc, params[ifree])
        if length(bestfit_val) != length(ifree)
            error("Length of best fit parameters ($(length(bestfit_val))) do not match number of free parameters ($(length(ifree)))")
        end

        pvalues[ifree] .= bestfit_val
        uncert[ifree] .= bestfit_unc
    end
    
    bestfit = OrderedDict{Symbol, FitComp}()
    ii = 0
    for (cname, comp) in model.comp
        fitcomp = OrderedDict{Symbol, Union{FitParam, Vector{FitParam}}}()
        accum = Vector{FitParam}()
        lastpname = :-
        for (pname, wpar) in getparams(comp)
            if (wpar.index == 0)  &&  (length(accum) > 0)
                fitcomp[lastpname] = deepcopy(accum)
                empty!(accum)
            end

            ii += 1
            params[ii].val = pvalues[ii] # Update model parameter values
            if wpar.index == 0
                fitcomp[pname] = FitParam(pvalues[ii], uncert[ii])
            else
                push!(accum, FitParam(pvalues[ii], uncert[ii]))
                lastpname = wpar.pname
            end
        end
        (length(accum) > 0)  &&  (fitcomp[lastpname] = deepcopy(accum))
        bestfit[cname] = FitComp(fitcomp)
    end

    result = FitResult(bestfit, length(data1d),
                       length(data1d) - length(ifree),
                       sum(abs2, (data1d.val .- model1d) ./ data1d.unc),
                       status, float(Base.time_ns() - elapsedTime) / 1.e9)
    return Wrap{FitResult}(result)
end


# ====================================================================
using LsqFit
mutable struct Minimizer <: AbstractMinimizer
end

support_param_limits(f::Minimizer) = false

function minimize(minimizer::Minimizer, evaluate::Function,
                  measure::Vector{Float64}, uncert::Vector{Float64},
                  params::Vector{Parameter})

    function callback(dummy::Vector{Float64}, pvalues::Vector{Float64})
        return evaluate(pvalues)
    end

    dom = collect(1.:length(measure))
    bestfit = LsqFit.curve_fit(callback, dom, measure, 1. ./ uncert, getfield.(params, :val))

    # Prepare output
    status = :NonOptimal
    if bestfit.converged
        status = :Optimal
    end

    error = LsqFit.margin_error(bestfit, 0.6827)
    return (status, getfield.(Ref(bestfit), :param), error)
end


