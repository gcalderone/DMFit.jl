module DataFitting

using Printf
using Statistics
using DataStructures

# ====================================================================
export test_component,
    Domain, CartesianDomain, getaxismin, getaxismax, getaxisextrema,
    Measures, FuncWrap, SimpleParam, flatten,
    Model, domain, evalcounter, resetcounters!,
    prepare!, evaluate!, fit

# ====================================================================
import Base.show
import Base.ndims
import Base.size
import Base.push!
import Base.length
import Base.getindex
import Base.propertynames
import Base.getproperty

include("Types.jl")


# ####################################################################
# Functions
#

# --------------------------------------------------------------------
"""
# Model

Constructor for the `Model` structure.
"""
Model(args::Vararg{Pair{Symbol, T}, N}) where {T<:AbstractComponent, N} = push!(Model(), args...)
function push!(model::Model, args::Vararg{Pair{Symbol, T}, N}) where {T<:AbstractComponent, N}
    for c in args
        getfield(model, :comp)[c[1]] = c[2]
    end
    return model
end

compcount(model::Model) = length(getfield(model, :comp))
components(model::Model) = getfield(model, :comp)
compiled(model::Model) = getfield(model, :compiled)

function compiled(model::Model, id::Int)
    c = compiled(model)
    @assert length(c) >= 1 "No model has been compiled"
    @assert 1 <= id <= length(c) "Invalid index (allowed range: 1 : " * string(length(c)) * ")"
    return c[id]
end

propertynames(model::Model) = collect(keys(getfield(model, :comp)))
getproperty(model::Model, s::Symbol) = get(getfield(model, :comp), s, nothing)

propertynames(res::BestFit) = collect(keys(getfield(res, :comp)))
getproperty(res::BestFit, s::Symbol) = get(getfield(res, :comp), s, nothing)

propertynames(res::BestFitComponent) = collect(keys(getfield(res, :params)))
getproperty(res::BestFitComponent, s::Symbol) = get(getfield(res, :params), s, nothing)



# --------------------------------------------------------------------
"""
# paramcount

Returns number of params in a component.
"""
paramcount(comp::AbstractComponent) = length(getparams(comp))

function getparams(comp::AbstractComponent)
    out = OrderedDict{Symbol, Parameter}()
    for pname in fieldnames(typeof(comp))
        tmp = getfield(comp, pname)
        if typeof(tmp) == Parameter
            out[pname] = tmp
        elseif typeof(tmp) == Vector{Parameter}
            for i in 1:length(tmp)
                out[Symbol(pname, i)] = tmp[i]
            end
        end
    end
    return out
end

function getparams(model::Model)
    out = OrderedDict{Symbol, Parameter}()
    for (cname, comp) in components(model)
        for (pname, par) in getparams(comp)
            out[Symbol(cname, "__", pname)] = par
        end
    end
    return out
end

getparamvalues(v::Union{Model, AbstractComponent}) = [par.val for (pname, par) in getparams(v)]

# --------------------------------------------------------------------
compdata(domain::AbstractDomain, comp::AbstractComponent) =
    error("Component " * string(typeof(comp)) * " must implement its own version of `compdata`.")


# --------------------------------------------------------------------
function CompiledExpression(model::Model, domain::AbstractDomain, exprs::Vector{Expr})
    function parse_model_expr(expr::Union{Symbol, Expr}, cnames, accum=Vector{Symbol}())
        if typeof(expr) == Expr
            # Parse the expression to check which components are involved
            for i in 1:length(expr.args)
                arg = expr.args[i]

                if typeof(arg) == Symbol
                    if arg in cnames # if it is one of the model components...
                        if i == 1  &&  expr.head == :call # ... and it is a function call...
                            error("Composite components are not yet supported")
                        else
                            push!(accum, arg)
                        end
                    end
                elseif typeof(arg) == Expr
                    parse_model_expr(arg, cnames, accum)
                else
                    println("Not handled: " * string(typeof(arg)))
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
        compinvolved = unique([compinvolved; parse_model_expr(expr, keys(components(model)))])
    end

    # TODO: Sort involved components according to insertion order
    kk = keys(components(model))
    sort!(compinvolved, lt=(a, b) -> findall(a .== kk)[1] < findall(b .== kk)[1])
    
    # Prepare the code for model evaluation
    code = Vector{String}()
    tmp = ""
    for (cname, comp) in components(model)
        for (pname, par) in getparams(comp)
            tmp *= ", $(cname)__$(pname)::Float64"
        end
    end
    push!(code, "(_ce::CompiledExpression $tmp, _unused_...) -> begin")
    # The last argument, _unused_, allows to push! further
    # components in the model after an expression has already been
    # prepared
    for (cname, comp) in components(model)
        if cname in compinvolved
            i = findall(cname .== compinvolved)[1]
            tmp = ""
            for (pname, par) in getparams(comp)
                (par.expr != "")  &&  (push!(code, "    $(cname)__$(pname) = " *
                                             replace(par.expr, "this__" => "$(cname)__")))
                tmp *= ", $(cname)__$(pname)"
            end
            push!(code, "  $cname = _evaluate!(_ce.compevals[$i], _ce.domain $tmp)")
        end
    end

    results = Vector{Vector{Float64}}()
    for i in 1:length(exprs)
        push!(code, "  @. _ce.results[$i] = (" * string(exprs[i]) * ")")
        push!(results, Vector{Float64}(undef, length(domain)))
    end
    push!(code, "  _ce.counter += 1")
    push!(code, "  return nothing")
    push!(code, "end")
    funct = eval(Meta.parse(join(code, "\n")))

    compevals = Vector{CompEvaluation}()
    for (cname, comp) in components(model)
        if cname in compinvolved
            tmp = CompEvaluation(0, compdata(domain, comp),
                                 Vector{Float64}(undef, paramcount(comp)),
                                 Vector{Float64}(undef, length(domain)))
            tmp.lastParams .= NaN
            push!(compevals, tmp)
        end
    end

    return CompiledExpression(join(code, "\n"), funct, 0,
                              deepcopy(domain), deepcopy(exprs), results,
                              compinvolved, compevals)
end
    
# --------------------------------------------------------------------
function _evaluate!(c::CompEvaluation, d::AbstractDomain, args...)
    if c.counter == 0
        c.counter += 1
        evaluate!(c.result, d, c.cdata, args...)
    else
        for i in 1:length(args)
            if c.lastParams[i] != args[i]
                c.lastParams .= args
                c.counter += 1
                evaluate!(c.result, d, c.cdata, args...)
            end
        end
    end
    return c.result
end

_evaluate!(ce::CompiledExpression, pvalues::Vector{Float64}) = Base.invokelatest(ce.funct, ce, pvalues...)

function evaluate!(model::Model, pvalues::Vector{Float64})
    for ce in compiled(model)
        _evaluate!(ce, pvalues)
    end
    return model
end

evaluate!(model::Model) = evaluate!(model, getparamvalues(model))


# --------------------------------------------------------------------
"""
# prepare!

Prepare a model to be evaluated on the given domain, with the given
mathematical expression.
"""
function prepare!(model::Model, domain::AbstractDomain, exprs::Vector{Expr})
    push!(compiled(model), CompiledExpression(model, domain, exprs))
    evaluate!(model)
    return model
end
prepare!(model::Model, domain::AbstractDomain, expr::Expr) = prepare!(model, domain, [expr])
prepare!(model::Model, domain::AbstractDomain, s::Symbol)  = prepare!(model, domain, [:(+$s)])


# --------------------------------------------------------------------
function getceval(ce::CompiledExpression, cname::Symbol)
    i = findall(ce.compnames .== cname)
    @assert length(i) == 1 "No component named $cname involved in compiled expression"
    return ce.compevals[i[1]]
end

"""
# domain

Return the domain associated to a model.
"""
domain(model::Model, id=1) = compiled(model, id).domain


"""
Returns a component evaluation.
"""
function getindex(model::Model, id::Int=1)
    out = compiled(model, id).results
    (length(out) == 1)  &&  (return out[1])
    return out
end
getindex(model::Model, id::Int, cname::Symbol) = getceval(compiled(model, id), cname).result

# --------------------------------------------------------------------
"""
# evalcounter

Return the number of times the model has been evaluated.
"""
evalcounter(model::Model, id::Int=1) = compiled(model, id).counter

"""
# evalcounter

Return the number of times a component has been evaluated.
"""
evalcounter(model::Model, id::Int, cname::Symbol) = getceval(compiled(model, id), cname).counter

"""
# resetcounters!

Reset model and components evaluation counters.
"""
function resetcounters!(model::Model)
    for ce in compiled(model)
        ce.counter = 0
        for ceval in ce.compevals
            ceval.counter = 0
        end
    end
end


# --------------------------------------------------------------------
function test_component(domain::AbstractLinearDomain, comp::AbstractComponent, iter=1)
    model = Model(:test => comp)
    prepare!(model, domain, :(+test))

    printstyled(color=:magenta, bold=true, "First evaluation:\n")
    @time result = evaluate!(model)
    if iter > 0
        println()
        printstyled(color=:magenta, bold=true, "Further evaluations ($iter):\n")

        @time begin
            for i in 1:iter
                result = evaluate!(model)
                compiled(model, 1).compevals[1].lastParams[1] = NaN  # Force re-calculation
            end
        end
    end
    show(model)
    return nothing
end
test_component(domain::AbstractCartesianDomain, comp::AbstractComponent, iter=1) =
    test_component(flatten(domain), comp, iter)


# ####################################################################
# Minimizer
#
support_param_limits(f::AbstractMinimizer) = false

"""
# fit

Fit a model against data, using the specified minimizer.
"""
function fit(model::Model, data::Vector{T}; minimizer=Minimizer()) where T<:AbstractMeasures
    elapsedTime = Base.time_ns()

    @assert typeof(minimizer) <: AbstractMinimizer
    @assert length(compiled(model)) >= 1

    # Check if the minimizer supports bounded parameters
    dict_params = getparams(model)
    pnames = collect(keys(dict_params))
    params = collect(values(dict_params))
    
    pvalues = getfield.(params, :val)
    ifree = findall(.! getfield.(params, :fixed))
    @assert length(ifree) > 0 "No free parameter in the model"

    if !support_param_limits(minimizer)
        if  (length(findall(isfinite.(getfield.(params, :low )))) > 0)  ||
            (length(findall(isfinite.(getfield.(params, :high)))) > 0)
            @warn "Parameter bounds are not supported by " * string(typeof(minimizer))
        end
    end

    # Prepare 1D arrays containing all the data and model results
    c1d_measure = Vector{Float64}()
    c1d_uncert  = Vector{Float64}()
    c1d_len  = Vector{Int}()
    push!(c1d_len, 0)
    for d in data
        append!(c1d_measure, d.measure)
        append!(c1d_uncert , d.uncert)
        push!(c1d_len, c1d_len[end] + length(d))
    end
    c1d_results = fill(0., length(c1d_measure))

    # Inner function to evaluate all the models and store the result in a 1D array
    all_ce = compiled(model)
    function evaluate1D(freepvalues::Vector{Float64})
        pvalues[ifree] .= freepvalues
        evaluate!(model, pvalues)

        ii = 1
        for ce in all_ce
            for v in ce.results
                c1d_results[c1d_len[ii]+1:c1d_len[ii+1]] .= v
                ii += 1
            end
        end
        return c1d_results
    end

    (status, bestfit_val, bestfit_unc) = minimize(minimizer, evaluate1D, c1d_measure, c1d_uncert, params[ifree])
    if length(bestfit_val) != length(ifree)
        error("Length of best fit parameters ($(length(bestfit_val))) do not match number of free parameters ($(length(ifree)))")
    end

    pvalues[ifree] .= bestfit_val
    uncert = fill(NaN, length(pvalues))
    uncert[ifree] .= bestfit_unc

    bestfit = OrderedDict{Symbol, BestFitComponent}()
    ii = 1
    for (cname, comp) in components(model)
        tmp = OrderedDict{Symbol, BestFitParameter}()
        for (pname, par) in getparams(comp)
            tmp[pname] = BestFitParameter(pvalues[ii], uncert[ii])
            ii += 1
        end
        bestfit[cname] = BestFitComponent(tmp)
    end
    result = FitResult(deepcopy(minimizer), BestFit(bestfit),
                       length(c1d_measure),
                       length(c1d_measure) - length(ifree),
                       sum(abs2, (c1d_measure .- c1d_results) ./ c1d_uncert),
                       status, float(Base.time_ns() - elapsedTime) / 1.e9)
    return result
end
fit(model::Model, data::AbstractData; minimizer=Minimizer()) = fit(model, [data]; minimizer=minimizer)


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
