# ====================================================================
#                         PRIVATE FUNCTIONS
# ====================================================================

# ____________________________________________________________________
# Model
#
function getparams(model::Model)
    out = OrderedDict{Symbol, Parameter}()
    for (cname, wcomp) in model.comp
        for (pname, par) in getparams(wcomp)
            out[Symbol(cname, compsep, pname)] = par
        end
    end
    return out
end

getparamvalues(v::Union{Model, AbstractComponent}) = [par.val for (pname, par) in getparams(v)]
function setparamvalues!(v::Union{Model, AbstractComponent}, pval::Vector{Float64})
    i = 0
    for (pname, par) in getparams(v)
        i += 1
        par.val = pval[i]
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
    model.buffer1d = Vector{Float64}(undef, model.index1d[end])
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
            push!(out, d1)
            ii += 1
        end
    end
    return out
end


# ____________________________________________________________________
# Component
#
function getparams(wcomp::WComponent)
    out = OrderedDict{Symbol, Parameter}()
    for pname in fieldnames(typeof(wcomp.comp))
        par = getfield(wcomp.comp, pname)
        if typeof(par) == Parameter
            par.cfixed = (par.fixed  ||  wcomp.fixed)
            if par.log
                par.val  = log10(par.val)
                par.low  = log10(par.low)
                par.high = log10(par.high)
            end
            out[pname] = par
        elseif typeof(par) == Vector{Parameter}
            for i in 1:length(par)
                par[i].cfixed = (par[i].fixed  ||  wcomp.fixed)
                out[Symbol(pname, i)] = par[i]
            end
        end
    end
    return out
end


# ____________________________________________________________________
# Parameter
#
isequal(a::Parameter, b::Parameter) = ((a._private.cname == b._private.cname)  &&
                                       (a._private.pname == b._private.pname)  &&
                                       (a._private.index == b._private.index))


# ____________________________________________________________________
# Instruments
#
# Instrument constructors
Instrument(dom::AbstractDomain) =
    Instrument(false, flatten(dom), dom, "", ()->nothing, 0,
               Vector{Symbol}(), Vector{CompEvaluation}(),
               Vector{Symbol}(), Vector{Expr}(), Vector{Bool}(), Vector{Vector{Float64}}())

function Instrument(model::Model, domain::AbstractDomain,
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
    for (cname, wcomp) in model.comp
        for (pname, par) in getparams(wcomp)
            tmp *= ", $(cname)$(compsep)$(pname)::Float64"
        end
    end
    push!(code, "(_instr::Instrument $tmp, _unused_...) -> begin")
    push!(code, "  domain = _instr.ldomain")
    # The last argument, _unused_, allows to push! further
    # components in the model after an expression has already been
    # prepared
    for (cname, wcomp) in model.comp
        if cname in compinvolved
            i = findall(cname .== compinvolved)[1]
            tmp = ""
            for (pname, par) in getparams(wcomp)
                (par.expr != "")  &&  (push!(code, "  $(cname)$(compsep)$(pname) = " *
                                             replace(par.expr, "this$(compsep)" => "$(cname)$(compsep)")))
                tmp *= ", $(cname)$(compsep)$(pname)"
            end
            if isfunction(wcomp.comp)
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

    compevals = CompEvals(model, domain, compinvolved)
    
    instr = Instrument(false, flatten(domain), domain, join(code, "\n"), funct, 0,
                       compinvolved, compevals,
                       deepcopy(labels), deepcopy(exprs), cmp, exprevals)
    return instr
end

function CompEvals(model::Model, domain::AbstractDomain, compinvolved::Vector{Symbol})
    compevals = Vector{CompEvaluation}()
    for (cname, wcomp) in model.comp
        if cname in compinvolved
            cd = cdata(wcomp.comp, domain)
            npar = length(getparams(wcomp))
            tmp = CompEvaluation(true, npar, fill(false, npar), 0, cd,
                                 Vector{Float64}(undef, npar),
                                 Vector{Float64}(undef, outsize(cd, domain)), [NaN])
            tmp.lastParams .= NaN
            push!(compevals, tmp)
        end
    end
    return compevals
end


# ____________________________________________________________________
# Expressions
#
newexprlabel(model::Model, id::Int) = Symbol(:expr, length(model.instruments[id].exprs)+1)
newexprlabel(model::Model, id::Int, n::Int) = Symbol.(Ref(:expr), length(model.instruments[id].exprs) .+ collect(1:n))

addexpr!(model::Model                                 , expr::Expr         ; cmp=true) = addexpr!(model,  1, [newexprlabel(model, 1)]              , [expr]       ; cmp=[cmp])
addexpr!(model::Model                                 , symbol::Symbol     ; cmp=true) = addexpr!(model,  1, [newexprlabel(model, 1)]              , [:(+$symbol)]; cmp=[cmp])
addexpr!(model::Model                                 , exprs::Vector{Expr}; cmp=true) = addexpr!(model,  1,  newexprlabel(model, 1, length(exprs)), exprs        ; cmp=fill(cmp, length(exprs)))

addexpr!(model::Model, id::Int                        , expr::Expr         ; cmp=true) = addexpr!(model, id, [newexprlabel(model, 1)]              , [expr]       ; cmp=[cmp])
addexpr!(model::Model, id::Int                        , symbol::Symbol     ; cmp=true) = addexpr!(model, id, [newexprlabel(model, 1)]              , [:(+$symbol)]; cmp=[cmp])
addexpr!(model::Model, id::Int                        , exprs::Vector{Expr}; cmp=true) = addexpr!(model, id,  newexprlabel(model, 1, length(exprs)), exprs        ; cmp=fill(cmp, length(exprs)))

addexpr!(model::Model,          label::Symbol         , expr::Expr         ; cmp=true) = addexpr!(model,  1, [label]                               , [expr]       ; cmp=[cmp])
addexpr!(model::Model,          label::Symbol         , symbol::Symbol     ; cmp=true) = addexpr!(model,  1, [label]                               , [:(+$symbol)]; cmp=[cmp])
addexpr!(model::Model,          labels::Vector{Symbol}, exprs::Vector{Expr}; cmp=true) = addexpr!(model,  1,  labels                               , exprs        ; cmp=fill(cmp, length(exprs)))

function addexpr!(model::Model, id::Int, labels::Vector{Symbol}, exprs::Vector{Expr}; cmp=Vector{Bool}())
    instr = model.instruments[id]
    append!(instr.exprnames, labels)
    append!(instr.exprs, exprs)
    append!(instr.exprcmp, cmp)
    instr.compile = true
    return labels[end]
end


# ____________________________________________________________________
function _evaluate!(c::CompEvaluation, d::AbstractDomain, args...)
    if c.fixed  &&  (c.counter > 0)
        return c.result
        isfinite(c.value[1])  &&  (return c.value)
    else
        @assert c.npar == length(args)
        curParams = convert(Vector{Float64}, [args...])
        neweval = false
        for i in 1:c.npar
            (c.lastParams[i] != curParams[i])  &&  (neweval = true)
            (c.log[i])  &&  (curParams[i] = 10. ^curParams[i])
        end
        if neweval
            c.lastParams .= args
            c.counter += 1
            evaluate!(c.cdata, c.result, d, curParams...)
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
end


function _evaluate!(model::Model)
    for id in 1:length(model.instruments)
        instr = model.instruments[id]
        if instr.compile
            tmp = Instrument(model, instr.domain, instr.exprnames, instr.exprs, instr.exprcmp)
            deleteat!(model.instruments, id)
            insert!(model.instruments, id, tmp)
            _evaluate!(model)
            prepareindex1D!(model)
        end
    end

    for instr in model.instruments
        for ii in 1:length(instr.compnames)
            wcomp = model.comp[instr.compnames[ii]]
            instr.compevals[ii].fixed = wcomp.fixed
            jj = 1
            for (pname, par) in getparams(wcomp)
                instr.compevals[ii].log[jj] = par.log
            end
        end
    end

    _evaluate!(model, getparamvalues(model))

    for instr in model.instruments
        for ceval in instr.compevals
            ceval.value[1] = NaN
            if ceval.fixed
                mm = extrema(ceval.result)
                (mm[1] == mm[2])  &&  (ceval.value[1] = mm[1])
            end
        end
    end
    return nothing
end


# ____________________________________________________________________
# Minimizer
#
function pvalues2FitComp(model::Model, pvalues::Vector{Float64}, uncert::Vector{Float64})
    bestfit = OrderedDict{Symbol, FitComp}()
    ii = 0
    for (cname, comp) in model.comp
        fitcomp = OrderedDict{Symbol, Union{FitParam, Vector{FitParam}}}()
        accum = Vector{FitParam}()
        lastpname = :-
        for (pname, par) in getparams(comp)
            if (par._private.index == 0)  &&  (length(accum) > 0)
                fitcomp[lastpname] = deepcopy(accum)
                empty!(accum)
            end

            ii += 1
            par.val = pvalues[ii] # Update model parameter values
            par._private.fitval = pvalues[ii]
            par._private.fitunc = uncert[ii]
            if par._private.index == 0
                fitcomp[pname] = FitParam(pvalues[ii], uncert[ii])
            else
                push!(accum, FitParam(pvalues[ii], uncert[ii]))
                lastpname = par._private.pname
            end
        end
        (length(accum) > 0)  &&  (fitcomp[lastpname] = deepcopy(accum))
        bestfit[cname] = FitComp(fitcomp)
    end
    return bestfit
end


function residuals1d(model::Model, data1d::Vector{Measures_1D}) where T<:AbstractMeasures
    ii = 1
    for id in 1:length(model.instruments)
        for jj in 1:length(model.instruments[id].exprnames)
            (model.instruments[id].exprcmp[jj])  ||  (continue)
            i1 = model.index1d[ii]+1
            i2 = model.index1d[ii+1]
            model.buffer1d[i1:i2] .=
                ((model.instruments[id].exprevals[jj] .- data1d[ii].val) ./ data1d[ii].unc)
            ii += 1
        end
    end
    return model.buffer1d
end


function _fit(model::Model, data::Vector{T}; kw...) where T<:AbstractMeasures
    pval = getparamvalues(model)
    res = _fit!(model, data; kw...)
    setparamvalues!(model, pval)
    _evaluate!(model)
    return res
end


function _fit!(model::Model, data::Vector{T}; dry=false, minimizer=lsqfit()) where T<:AbstractMeasures
    elapsedTime = Base.time_ns()

    @assert typeof(minimizer) <: AbstractMinimizer
    @assert length(model.instruments) >= 1
    _evaluate!(model)
    params = collect(values(getparams(model)))
    pvalues = getfield.(params, :val)
    uncert = fill(NaN, length(pvalues))
    ifree = findall(.! getfield.(params, :cfixed))
    @assert length(ifree) > 0 "No free parameter in the model"

    # Inner function to evaluate all the models and store the result in a 1D array
    data1d = data1D(model, data)
    function eval_residuals1d(freepvalues::Vector{Float64},
                              model=model, pvalues=pvalues, ifree=ifree, data1d=data1d)
        pvalues[ifree] .= freepvalues
        _evaluate!(model, pvalues)
        return residuals1d(model, data1d)
    end
    eval_residuals1d(getfield.(params[ifree], :val))
    #Main.code_warntype(eval_residuals1d, (Vector{Float64},))

    status = :NonOptimal
    if !dry
        (status, bestfit_val, bestfit_unc) = minimize(minimizer, eval_residuals1d, params[ifree])
        @assert length(bestfit_val) == length(ifree)
        pvalues[ifree] .= bestfit_val
        uncert[ifree] .= bestfit_unc
        eval_residuals1d(pvalues[ifree]) # Re-evaluate model to ensure best fit values are used
    end

    cost = sum(abs2, model.buffer1d)
    dof = length(model.buffer1d) - length(ifree)
    result = FitResult(pvalues2FitComp(model, pvalues, uncert),
                       length(model.buffer1d), dof,
                       cost, status, logccdf(Chisq(dof), cost) * log10(exp(1)),
                       float(Base.time_ns() - elapsedTime) / 1.e9)
    return UI{FitResult}(result)
end


# ====================================================================
using LsqFit
mutable struct lsqfit <: AbstractMinimizer
end

function minimize(minimizer::lsqfit, func::Function, params::Vector{Parameter})
    ndata = length(func(getfield.(params, :val)))
    bestfit = LsqFit.curve_fit((dummy, pvalues) -> func(pvalues),
                               1.:ndata, fill(0., ndata),
                               getfield.(params, :val),
                               lower=getfield.(params, :low),
                               upper=getfield.(params, :high))
    status = :NonOptimal
    (bestfit.converged)  &&  (status = :Optimal)
    error = LsqFit.margin_error(bestfit, 0.6827)
    return (status, getfield.(Ref(bestfit), :param), error)
end



macro enable_CMPFit()
    return esc(:(
        import DataFitting.minimize;

        mutable struct cmpfit <: DataFitting.AbstractMinimizer;
        config::CMPFit.Config;
        cmpfit() = new(CMPFit.Config());
        end;

        function minimize(minimizer::cmpfit, func::Function, params::Vector{DataFitting.Parameter});
        guess = getfield.(params, :val);
        low   = getfield.(params, :low);
        high  = getfield.(params, :high);
        parinfo = CMPFit.Parinfo(length(guess));
        for i in 1:length(guess);
        llow  = isfinite(low[i])   ?  1  :  0;
        lhigh = isfinite(high[i])  ?  1  :  0;
        parinfo[i].limited = (llow, lhigh);
        parinfo[i].limits  = (low[i], high[i]);
        end;
        bestfit = CMPFit.cmpfit((pvalues) -> func(pvalues),
                                guess, parinfo=parinfo, config=minimizer.config);
        return (:Optimal, getfield.(Ref(bestfit), :param), getfield.(Ref(bestfit), :perror));
        end;
    ))
end
