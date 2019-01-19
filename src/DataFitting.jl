__precompile__(true)

module DataFitting

using Printf
using Statistics

# ====================================================================
export test_component,
    Domain, CartesianDomain, getaxismin, getaxismax, getaxisextrema,
    Measures, FuncWrap, SimpleParam,
    Model, prepare!, evaluate!, setcompvalue!, resetcounters!, fit!,
    domain, evalcounter

# ====================================================================
import Base.push!
import Base.show
import Base.ndims
import Base.size
import Base.reshape
import Base.length
import Base.iterate
import Base.keys
import Base.getindex
import Base.setindex!


# ====================================================================
include("HashVector.jl")
include("Types.jl")


# ####################################################################
# Functions
#

# --------------------------------------------------------------------
"""
# nparams

Returns number of params in a component.
"""
function nparams(comp::AbstractComponent)
    count = 0
    for pname in fieldnames(typeof(comp))
        if fieldtype(typeof(comp), pname) == Parameter
            count += 1
        end
        if fieldtype(typeof(comp), pname) == Vector{Parameter}
            count += length(getfield(comp, pname))
        end
    end
    return count
end


# --------------------------------------------------------------------
"""
# Model

Constructor for the `Model` structure.
"""
function Model(c::Pair{Symbol, T}, args...) where T<:AbstractComponent
    model = Model()
    push!(model.comp,      c[1], c[2])
    push!(model.altValues, c[1], NaN)

    for c in args
        @assert typeof(c) <: Pair
        @assert typeof(c[1]) == Symbol
        @assert typeof(c[2]) <: AbstractComponent
        push!(model.comp,      c[1], c[2])
        push!(model.altValues, c[1], NaN)
    end
    return model
end

# Forward methods to Model.comp
length(model::Model) = length(model.comp)
keys(model::Model) = keys(model.comp)
getindex(model::Model, key::Symbol) = getindex(model.comp, key)

iterate(model::Model) = iterate(model.comp)
iterate(model::Model, i::Int) = iterate(model.comp, i)

function push!(model::Model, cname::Symbol, comp::T) where T <: AbstractComponent
    push!(model.comp, cname, comp)
    push!(model.altValues, cname, NaN)
    return model
end


function getparams(comp::T, cname::String="") where T <: AbstractComponent
    prefix = cname
    (prefix != "")  &&  (prefix *= "__")
    pnames = Vector{String}()
    params = Vector{Parameter}()
    count = 0
    for pname in fieldnames(typeof(comp))
        isVector = false
        if fieldtype(typeof(comp), pname) == Parameter
            push!(params, getfield(comp, pname))
            wname = prefix * string(pname)
            push!(pnames, wname)
        elseif fieldtype(typeof(comp), pname) == Vector{Parameter}
            isVector = true
            count = 0
            push!(params, getfield(comp, pname)...)
            for j in 1:length(getfield(comp, pname))
                wname = prefix * string(pname) * string(j)
                push!(pnames, wname)
            end
        end
    end
    return (pnames, params)
end


function getsetparams!(model::Model, newvalues=Vector{Float64}(), parToUpdate=Vector{Symbol}())
    out_cnames = Vector{Symbol}()
    out_pnames = Vector{Symbol}()
    out_wnames = Vector{Symbol}()
    out_params = Vector{Parameter}()

    count = 0
    for (cname, comp) in model.comp
        (pnames, params) = getparams(comp)
        for i in 1:length(params)
            par = params[i]
            pname = pnames[i]
            wname = string(cname) * "__" * pname

            newval = NaN
            if length(newvalues) != 0
                if length(parToUpdate) != 0
                    j = findall(wname .== parToUpdate)
                    (length(j) == 1)  &&  (newval = newvalues[j[1]])
                else
                    count += 1
                    newval = newvalues[count]
                end
            end
            if isfinite(newval)
                @assert (par.low <= newval <= par.high) "Value of parameter $(wname) is outside limits: $(par.low) < $(newval) < $(par.high)"
                par.val = newval
            end

            @assert (par.low <= par.val <= par.high) "Value of parameter $(wname) is outside limits: $(par.low) < $(par.val) < $(par.high)"
            (isfinite(model.altValues[cname]))  &&  (par.fixed = true)

            push!(out_cnames, cname)
            push!(out_pnames, Symbol(pname))
            push!(out_wnames, Symbol(wname))
            push!(out_params, par)
        end
    end
    return (out_cnames, out_pnames, out_wnames, out_params)
end



# --------------------------------------------------------------------
compdata(domain::AbstractDomain, comp::AbstractComponent) =
    error("Component " * string(typeof(comp)) * " must implement its own version of `compdata`.")



# --------------------------------------------------------------------
"""
# prepare!

Recompile the expressions after one (or more) parameter expressions have been changed
"""
function prepare!(model::Model)
    bkp = deepcopy(model)

    # Delete all compiled expressions and results
    empty = Model()
    model.compiled = empty.compiled
    model.results = empty.results
    model.domainid = empty.domainid

    # Recompile expressions
    for i in 1:length(bkp.compiled)
        prepare!(model, bkp.compiled[1].domain, bkp.compiled[1].exprs)
    end

    return model
end


# --------------------------------------------------------------------
"""
# prepare!

Prepare a model to be evaluated on the given domain, with the given
mathematical expression to join the components.
"""
function prepare!(model::Model, domain::AbstractDomain, exprs::Vector{Expr}; cache=true)

    function CompiledExpression1(model::Model, exprs::Vector{Expr}; cache=true)
        function parse_model_expr(expr, cnames, simple=Vector{Symbol}(), composite=Vector{Symbol}())
            if typeof(expr) == Expr
                # Parse the expression to check which components are involved
                for i in 1:length(expr.args)
                    arg = expr.args[i]
                    
                    if typeof(arg) == Symbol
                        if arg in cnames # if it is one of the model components...
                            if i == 1  &&  expr.head == :call # ... and it is a function call...
                                push!(composite, arg)
                            else
                                push!(simple, arg)
                            end
                        end
                    elseif typeof(arg) == Expr
                        parse_model_expr(arg, cnames, simple, composite)
                    else
                        println("Not handled: " * string(typeof(arg)))
                    end
                end
            else
                push!(simple, expr)
            end
            return (simple, composite)
        end

        # Check which components are involved.  TODO: handle composite components
        compInvolved = Vector{Symbol}()
        for i in 1:length(exprs)
            (a, b) = parse_model_expr(exprs[i], keys(model))
            push!(compInvolved, a...)
        end
        compInvolved = unique(compInvolved)
        
        # Prepare the code for model evaluation
        code = Vector{String}()
        (cnames, pnames, wnames, params) = getsetparams!(model)
        push!(code, "(_m::CompiledExpression, _altValues::HashVector{Float64}, _results::Vector{Vector{Float64}}, " *
              join(string.(wnames) .* "::Float64", ", ") * ", _unused_...) -> begin")
        # The last argument, _unused_, allows to push! further
        # components in the model after an expression has already been
        # prepared
        countAll = 0
        countInv = 0
        for (cname, comp) in model
            countAll += 1
            (!(cname in compInvolved))  &&  (continue)
            countInv += 1

            push!(code, "  $cname = _m.cevals[$countInv].result")
            push!(code, "  if isfinite(_altValues[$countAll])")
            push!(code, "    $cname = _altValues[$countAll]")
            push!(code, "  else")

            (wnames, params) = getparams(comp, string(cname))
            for i in 1:length(params)
                par = params[i]
                wname = wnames[i]

                if par.expr != ""
                    push!(code, "    $(wname) = " * replace(par.expr, "this__" => "$(cname)__"))
                end
            end

            tmp2 = "      cached_evaluate!(_m.cevals[$countInv], _m.ldomain"
            for i in 1:length(params)
                wname = wnames[i]
                tmp2 *= ", $(wname)"
            end
            tmp2 *= ")"
            push!(code, tmp2)
            push!(code, "  end")
            push!(code, "")
        end
        
        j = 1 + length(model.results)
        push!(code, "  if length(_results[$j]) <= 1")
        for i in 1:length(exprs)
            j = i + length(model.results)
            push!(code, "    _results[$j] = @. (" * string(exprs[i]) * ")")
        end
        push!(code, "  else")
        for i in 1:length(exprs)
            j = i + length(model.results)
            push!(code, "    @.( _results[$j] = " * string(exprs[i]) * ")")
        end
        push!(code, "  end")
        push!(code, "")
        push!(code, "  _m.counter += 1")
        push!(code, "  return nothing")
        push!(code, "end")
        code = join(code, "\n")

        # Evaluate the code
        expr2 = Meta.parse(code)
        funct = eval(expr2)
        return (code, funct, compInvolved)
    end
    

    function CompiledExpression2(model::Model, domain::AbstractDomain, compInvolved::Vector{Symbol})
        ldomain = domain
        if typeof(domain) <: AbstractCartesianDomain
            ldomain = flatten(domain)  # <- TYPE INSTABILITY
        end

        cevals = HashVector{ComponentEvaluation}()
        for (cname, comp) in model
            (cname in compInvolved)  ||   (continue)
            buffer = Vector{Float64}(undef, length(ldomain))
            cdata = compdata(ldomain, comp)

            tmp = ComponentEvaluation(cdata, Vector{Float64}(undef, nparams(comp)), buffer, 0)
            tmp.lastParams .= NaN
            push!(cevals, cname, tmp)
        end

        return (ldomain, cevals)
    end

    (code, funct, compInvolved) = CompiledExpression1(model, exprs, cache=cache)
    (ldomain, cevals) = CompiledExpression2(model, domain, compInvolved)
    push!(model.compiled, CompiledExpression(code, funct, deepcopy(exprs), compInvolved, 0,
                                             ldomain, deepcopy(domain), cevals))
    for i in 1:length(exprs)
        push!(model.results, Vector{Float64}())
        push!(model.domainid, length(model.compiled))
    end
    evaluate!(model)
    return model
end
prepare!(model::Model, domain::AbstractDomain, expr::Expr; cache=true) = prepare!(model, domain, [expr], cache=cache)
prepare!(model::Model, domain::AbstractDomain, s::Symbol; cache=true) = prepare!(model, domain, [:(+$s)], cache=cache)



# --------------------------------------------------------------------
function checkID(model::Model, id::Int)
    @assert length(model.compiled) >= 1 "No model has been compiled"
    @assert 1 <= id <= length(model.compiled) "Invalid index (allowed range: 1 : " * string(length(model.compiled)) * ")"
end

"""
# evaluate!

Evaluate the model(s)
"""
function evaluate!(model::Model, id::Int, pvalues::Vector{Float64})
    checkID(model, id)
    return Base.invokelatest(model.compiled[id].funct, model.compiled[id], model.altValues, model.results, pvalues...)
end

function evaluate!(model::Model, pvalues::Vector{Float64})
    for id in 1:length(model.compiled)
        evaluate!(model, id, pvalues)
    end
    return model
end

function evaluate!(model::Model)
    (cnames, pnames, wnames, params) = getsetparams!(model)
    pvalues = getfield.(params, :val)
    return evaluate!(model, pvalues)
end


function cached_evaluate!(c::ComponentEvaluation, d::AbstractDomain, args...)
    if c.counter == 0
            c.counter += 1
            return evaluate!(c.result, d, c.cdata, args...)
    end
    for i in 1:length(args)
        if c.lastParams[i] != args[i]
            c.lastParams .= args
            c.counter += 1
            return evaluate!(c.result, d, c.cdata, args...)
        end
    end
    return c.result
end


# --------------------------------------------------------------------
"""
# setcompvalue!

Set a component value in the model
"""
function setcompvalue!(model::Model, cname::Symbol, value::Number=NaN)
    model.altValues[cname] = float(value)
    evaluate!(model)
    return value
end

"""
# domain

Return the domain associated to a model.
"""
function domain(model::Model, id=1)
    checkID(model, id)
    return model.compiled[id].domain
end



"""
Returns a component evaluation.
"""
function (model::Model)(comp::Symbol, id=1)
    checkID(model, id)

    if comp in keys(model.comp)
        if isfinite(model.altValues[comp])
            return model.altValues[comp]
        end

        if comp in model.compiled[id].compInvolved
            return model.compiled[id].cevals[comp].result
        end
        return [NaN]
    end
    error("No component named $comp is defined.")
end


"""
Returns a model evaluation with the parameter values given as input.
"""
function (model::Model)(id::Int, params::Vararg{Pair{Symbol,T}}) where T <: Number
    @assert 1 <= id <= length(model.results) "Invalid index (allowed range: 1 : " * string(length(model.results)) * ")"
    if length(params) > 0
        pnames = Vector{Symbol}()
        pvalues = Vector{Float64}()
        for (pname, pvalue) in params
            push!(pnames, pname)
            push!(pvalues, pvalue)
        end
        getsetparams!(model, pvalues, pnames)
        evaluate!(model)
    end
    return model(id)
end

(model::Model)() = model(1)
function (model::Model)(id::Int)
    @assert 1 <= id <= length(model.results) "Invalid index (allowed range: 1 : " * string(length(model.results)) * ")"
    return model.results[id]
end

function (model::Model)(params::Vararg{Pair{Symbol,T}}) where T <: Number
    return model(1, params...)
end

"""
# evalcounter

Return the number of times the model has been evaluated.
"""
function evalcounter(model::Model, id=1)
    checkID(model, id)
    return model.compiled[id].counter
end

"""
# evalcounter

Return the number of times a component has been evaluated.
"""
function evalcounter(model::Model, comp::Symbol, id=1)
    checkID(model, id)
    return model.compiled[id].cevals[comp].counter
end

"""
# resetcounters!

Reset model and components evaluation counters.
"""
function resetcounters!(model::Model)
    for id in 1:length(model.compiled)
        for j in 1:length(model.compiled[id].cevals)
            model.compiled[id].cevals[j].counter = 0
        end
        model.compiled[id].counter = 0
    end
end


# --------------------------------------------------------------------
function test_component(domain::AbstractLinearDomain, comp::AbstractComponent, iter=1)
    #println()
    #printstyled(color=:magenta, bold=true, "================================================================================\n")
    printstyled(color=:magenta, bold=true, "Profiling component:\n")
    model = Model(:test => comp)
    show(model)

    prepare!(model, domain, :(+test), cache=false)

    println()
    printstyled(color=:magenta, bold=true, "First evaluation:\n")
    @time result = evaluate!(model)

    if iter > 0
        println()
        printstyled(color=:magenta, bold=true, "Further evaluations ($iter):\n")

        @time begin
            for i in 1:iter
                result = evaluate!(model)
            end
        end
    end
    #println()
    #show(model)

    return nothing
end
test_component(domain::AbstractCartesianDomain, comp::AbstractComponent, iter=1) =
    test_component(flatten(domain), comp, iter)



# ####################################################################
# Minimizer
#
support_param_limits(f::AbstractMinimizer) = false

"""
# fit!

Fit a model against data, using the specified minimizer.
"""
function fit!(model::Model, data::Vector{T}; minimizer=Minimizer()) where T<:AbstractMeasures
    @assert typeof(minimizer) <: AbstractMinimizer
    elapsedTime = Base.time_ns()

    @assert length(model.comp) >= 1
    @assert length(model.compiled) >= 1
    @assert length(model.results) == length(data) "Model has " * string(length(model.results)) * " expressions but " * string(length(data)) * " data were given."

    # Check if the minimizer supports bounded parameters
    (cnames, pnames, wnames, params) = getsetparams!(model)
    pvalues = getfield.(params, :val)
    ifree = findall(.! getfield.(params, :fixed))
    @assert length(ifree) > 0 "No free parameter in the model"

    if !support_param_limits(minimizer)
        if  (length(findall(isfinite.(getfield.(params, :low )))) > 0)  ||
            (length(findall(isfinite.(getfield.(params, :high)))) > 0)
            printstyled(color=:red, bold=true, "Parameter bounds are not supported by " * string(typeof(minimizer)) * "\n")
        end
    end

    # Prepare 1D arrays containing all the data and model results
    c1d_measure = Vector{Float64}()
    c1d_uncert  = Vector{Float64}()
    c1d_len  = Vector{Int}()
    push!(c1d_len, 0)
    for i in 1:length(data)
        m = model.compiled[1]
        if length(model.compiled) > 1
            m = model.compiled[model.domainid[i]]
        end
        tmp = data[i]
        if ndims(tmp) > 1
            tmp = flatten(data[i], m.domain)
        end
        push!(c1d_measure, tmp.measure...)
        push!(c1d_uncert , tmp.uncert...)
        push!(c1d_len , c1d_len[end] + length(tmp))
    end
    c1d_results = fill(0., length(c1d_measure))

    # Inner function to evaluate all the models and store the result in a 1D array
    function evaluate1D(freepvalues::Vector{Float64})
        pvalues[ifree] .= freepvalues
        evaluate!(model, pvalues)

        for i in 1:length(data)
            rr = model.results[1]
            if length(model.results) > 1
                rr = model.results[i]
            end
            c1d_results[c1d_len[i]+1:c1d_len[i+1]] .= rr
        end
        return c1d_results
    end

    #try
    (status, bestfit_val, bestfit_unc) = minimize(minimizer, evaluate1D, c1d_measure, c1d_uncert, params[ifree])
    if length(bestfit_val) != length(ifree)
        error("Length of best fit parameters ($(length(bestfit_val))) do not match number of free parameters ($(length(ifree)))")
    end

    pvalues[ifree] .= bestfit_val
    uncert = fill(NaN, length(pvalues))
    uncert[ifree] .= bestfit_unc

    bestfit = HashVector{FitParameter}()
    for i in 1:length(params)
        push!(bestfit, wnames[i], FitParameter(pvalues[i], uncert[i]))
    end

    getsetparams!(model, pvalues)
    result = FitResult(deepcopy(minimizer), bestfit,
                       length(c1d_measure),
                       length(c1d_measure) - length(ifree),
                       sum(abs2, (c1d_measure .- c1d_results) ./ c1d_uncert),
                       status, float(Base.time_ns() - elapsedTime) / 1.e9)
    return result
    # catch err
    #     printstyled(color=:red, err)
    #     println()
    #     println("Call stack:")
    #     for s in stacktrace(catch_backtrace())
    #         println(s)
    #     end
    #     println()
    # end
end
fit!(model::Model, data::AbstractData; minimizer=Minimizer()) = fit!(model, [data]; minimizer=minimizer)



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

end
