module DataFitting

using Printf
using Statistics
using DataStructures

# ====================================================================
export test_component,
    Domain, CartesianDomain, getaxismin, getaxismax, getaxisextrema,
    Measures, FuncWrap, SimpleParam, flatten,
    Model, addcomponent!, addinstrument!, addexpr!, domain, evalcounter, resetcounters!,
    addinstrument!, evaluate!, fit

# ====================================================================
import Base.show
import Base.ndims
import Base.size
import Base.length
import Base.getindex
import Base.propertynames
import Base.getproperty

include("Types.jl")
include("show.jl")

const compsep = "_"
const showprefix = "    "

# ####################################################################
# Functions
#

# ____________________________________________________________________
"""
# Model

Constructor for the `Model` structure.
"""
Model(args::Vararg{Pair{Symbol, T}, N}) where {T<:AbstractComponent, N} = addcomponent!(Model(), args...)
function addcomponent!(model::Model, args::Vararg{Pair{Symbol, T}, N}) where {T<:AbstractComponent, N}
    for c in args
        getfield(model, :comp)[c[1]] = c[2]
    end
    return model
end

compcount(model::Model) = length(getfield(model, :comp))
components(model::Model) = getfield(model, :comp)
instruments(model::Model) = getfield(model, :instruments)

function instruments(model::Model, id::Int)
    c = instruments(model)
    @assert length(c) >= 1 "No model has been compiled"
    @assert 1 <= id <= length(c) "Invalid index (allowed range: 1 : " * string(length(c)) * ")"
    return c[id]
end

propertynames(model::Model) = collect(keys(getfield(model, :comp)))
getproperty(model::Model, s::Symbol) = get(getfield(model, :comp), s, nothing)

propertynames(res::BestFit) = collect(keys(getfield(res, :comp)))
getproperty(res::BestFit, s::Symbol) = get(getfield(res, :comp), s, nothing)

propertynames(res::BestFitComp) = collect(keys(getfield(res, :params)))
getproperty(res::BestFitComp, s::Symbol) = get(getfield(res, :params), s, nothing)



# ____________________________________________________________________
"""
# paramcount

Returns number of params in a component.
"""
paramcount(comp::AbstractComponent) = length(getparams(comp))

function getparams(comp::AbstractComponent)
    out = OrderedDict{Symbol, WrapParameter}()
    for pname in fieldnames(typeof(comp))
        par = getfield(comp, pname)
        if typeof(par) == Parameter
            out[pname] = WrapParameter(pname, 0, par)
        elseif typeof(par) == Vector{Parameter}
            for i in 1:length(par)
                out[Symbol(pname, i)] = WrapParameter(pname, i, par[i])
            end
        end
    end
    return out
end


function getparams(model::Model)
    out = OrderedDict{Symbol, WrapParameter}()
    for (cname, comp) in components(model)
        for (pname, par) in getparams(comp)
            out[Symbol(cname, compsep, pname)] = par
        end
    end
    return out
end

getparamvalues(v::Union{Model, AbstractComponent}) = [wpar.par.val for (pname, wpar) in getparams(v)]

# ____________________________________________________________________
compdata(domain::AbstractDomain, comp::AbstractComponent) =
    error("Component " * string(typeof(comp)) * " must implement its own version of `compdata`.")

# ____________________________________________________________________
newexprlabel(model::Model, id::Int) = Symbol(:expr, length(instruments(model, id).exprs)+1)
newexprlabel(model::Model, id::Int, n::Int) = Symbol.(Ref(:expr), length(instruments(model, id).exprs) .+ collect(1:n))

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
    (length(cmp) == 0)  && (cmp = fill(true, length(exprs)))
    @assert length(labels) == length(exprs)
    @assert length(cmp) == length(exprs)

    instrument = instruments(model, id)
    labels = [instrument.exprnames; labels]
    exprs = [instrument.exprs; exprs]
    cmp = [instrument.exprcmp; cmp]
    deleteat!(instruments(model), id)
    instrument = Instrument(instrument.label, instrument.domain)
    
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

    # Sort involved components according to insertion order
    kk = keys(components(model))
    sort!(compinvolved, lt=(a, b) -> findall(a .== kk)[1] < findall(b .== kk)[1])
    
    # Prepare the code for model evaluation
    code = Vector{String}()
    tmp = ""
    for (cname, comp) in components(model)
        for (pname, wpar) in getparams(comp)
            tmp *= ", $(cname)$(compsep)$(pname)::Float64"
        end
    end
    push!(code, "(_instrument::Instrument $tmp, _unused_...) -> begin")
    # The last argument, _unused_, allows to push! further
    # components in the model after an expression has already been
    # prepared
    for (cname, comp) in components(model)
        if cname in compinvolved
            i = findall(cname .== compinvolved)[1]
            tmp = ""
            for (pname, wpar) in getparams(comp)
                par = wpar.par
                (par.expr != "")  &&  (push!(code, "  $(cname)$(compsep)$(pname) = " *
                                             replace(par.expr, "this$(compsep)" => "$(cname)$(compsep)")))
                tmp *= ", $(cname)$(compsep)$(pname)"
            end
            push!(code, "  $cname = _evaluate!(_instrument.compevals[$i], _instrument.domain $tmp)")
        end
    end

    exprevals = Vector{Vector{Float64}}()
    for i in 1:length(exprs)
        push!(code, "  @. _instrument.exprevals[$i] = (" * string(exprs[i]) * ")")
        push!(exprevals, Vector{Float64}(undef, length(instrument.domain)))
    end
    push!(code, "  _instrument.counter += 1")
    push!(code, "  return nothing")
    push!(code, "end")
    funct = eval(Meta.parse(join(code, "\n")))

    compevals = Vector{CompEvaluation}()
    for (cname, comp) in components(model)
        if cname in compinvolved
            tmp = CompEvaluation(0, compdata(instrument.domain, comp),
                                 Vector{Float64}(undef, paramcount(comp)),
                                 Vector{Float64}(undef, length(instrument.domain)))
            tmp.lastParams .= NaN
            push!(compevals, tmp)
        end
    end

    instrument.code = join(code, "\n")
    instrument.funct = funct
    instrument.counter = 0
    instrument.compnames = compinvolved
    instrument.compevals = compevals
    instrument.exprnames = deepcopy(labels)
    instrument.exprs = deepcopy(exprs)
    instrument.exprcmp = cmp
    instrument.exprevals = exprevals
    insert!(instruments(model), id, instrument)
    evaluate!(model)
    return model
end
    
# ____________________________________________________________________
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

_evaluate!(instrument::Instrument, pvalues::Vector{Float64}) = Base.invokelatest(instrument.funct, instrument, pvalues...)

function evaluate!(model::Model, pvalues::Vector{Float64})
    for instrument in instruments(model)
        _evaluate!(instrument, pvalues)
    end
    return model
end

evaluate!(model::Model) = evaluate!(model, getparamvalues(model))


# ____________________________________________________________________
"""
# addinstrument!

Prepare a model to be evaluated on the given domain, with the given
mathematical expression.
"""

function addinstrument!(model::Model, domain::AbstractDomain; label="none")
    push!(instruments(model), Instrument(label, domain))
    return length(instruments(model))
end


# function addinstrument!(model::Model)
#     bkg = deepcopy(model)
#     empty!(getfield(model, :instruments))
#     for instrument in instruments(bkg)
#         addinstrument!(model, instrument.domain, instrument.exprs)
#     end
#     return model
# end


# ____________________________________________________________________
function getceval(instrument::Instrument, cname::Symbol)
    i = findall(instrument.compnames .== cname)
    @assert length(i) == 1 "No component named $cname involved in compiled expression"
    return instrument.compevals[i[1]]
end

"""
# domain

Return the domain associated to a model.
"""
domain(model::Model, id=1) = instruments(model, id).domain


"""
Returns a component evaluation.
"""
function getindex(model::Model, id::Int=1, expr::Int=1)
    out = instruments(model, id).exprevals
    @assert expr <= length(out) "Invalid expression index: $expr"
    return out[expr]
end
getindex(model::Model, id::Int, cname::Symbol) =
    getceval(instruments(model, id), cname).result

# ____________________________________________________________________
"""
# evalcounter

Return the number of times the model has been evaluated.
"""
evalcounter(model::Model, id::Int=1) = instruments(model, id).counter

"""
# evalcounter

Return the number of times a component has been evaluated.
"""
evalcounter(model::Model, id::Int, cname::Symbol) = getceval(instruments(model, id), cname).counter

"""
# resetcounters!

Reset model and components evaluation counters.
"""
function resetcounters!(model::Model)
    for instrument in instruments(model)
        instrument.counter = 0
        for ceval in instrument.compevals
            ceval.counter = 0
        end
    end
end


# ____________________________________________________________________
function test_component(domain::AbstractLinearDomain, comp::AbstractComponent, iter=1)
    model = Model(:test => comp)
    addinstrument!(model, domain, :(+test))

    printstyled(color=:magenta, bold=true, "First evaluation:\n")
    @time result = evaluate!(model)
    if iter > 0
        println()
        printstyled(color=:magenta, bold=true, "Further evaluations ($iter):\n")

        @time begin
            for i in 1:iter
                result = evaluate!(model)
                instruments(model, 1).compevals[1].lastParams[1] = NaN  # Force re-calculation
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
    @assert length(instruments(model)) >= 1

    # Check if the minimizer supports bounded parameters
    params = getfield.(values(getparams(model)), :par)
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
    all_instrument = instruments(model)
    function evaluate1D(freepvalues::Vector{Float64})
        pvalues[ifree] .= freepvalues
        evaluate!(model, pvalues)

        ii = 1
        for instrument in all_instrument
            for v in instrument.exprevals
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

    bestfit = OrderedDict{Symbol, BestFitComp}()
    ii = 1
    for (cname, comp) in components(model)
        bestfitcomp = OrderedDict{Symbol, Union{BestFitParam, Vector{BestFitParam}}}()
        accum = Vector{BestFitParam}()
        lastpname = :-
        for (pname, wpar) in getparams(comp)
            if (wpar.index == 0)  &&  (length(accum) > 0)
                bestfitcomp[lastpname] = deepcopy(accum)
                empty!(accum)
            end

            tmp = BestFitParam(pvalues[ii], uncert[ii]); ii += 1
            if wpar.index == 0
                bestfitcomp[pname] = tmp
            else
                push!(accum, tmp)
                lastpname = wpar.pname
            end
        end
        if length(accum) > 0
            bestfitcomp[lastpname] = deepcopy(accum)
        end
        bestfit[cname] = BestFitComp(bestfitcomp)
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
