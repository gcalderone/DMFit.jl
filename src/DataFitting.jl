module DataFitting

using Printf, PrettyTables
using Statistics, Distributions
using DataStructures
using LsqFit

import Base.push!
import Base.show
import Base.ndims
import Base.size
import Base.length
import Base.getindex
import Base.reshape


export Domain, CartesianDomain, Measures,
    Prediction, Model, evaluate, parindex, thaw, freeze, fit!


include("domain.jl")


# ====================================================================
# Parameter
#
mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    step::Float64
    free::Bool
    Parameter(value::Number) = new(float(value), -Inf, +Inf, NaN, true)
end

# ====================================================================
# A *component* is a generic implementation of a constituent part of a
# model.
#
# A component must inherit `AbstractComponent`, and implement the
# `ceval_data` and `evaluate` methods
abstract type AbstractComponent end


function getparams(comp::AbstractComponent)
    params = OrderedDict{Symbol, Parameter}()
    for pname in fieldnames(typeof(comp))
        par = getfield(comp, pname)
        if isa(par, Parameter)
            params[pname] = par
        elseif isa(par, Vector{Parameter})
            for i in 1:length(par)
                params[Symbol(pname, i)] = par[i]
            end
        end
    end
    params
end

# ====================================================================
# CompEval: a wrapper for a component evaluated on a specific domain
#
mutable struct CompEval{TDomain <: AbstractDomain, TComp <: AbstractComponent}
    domain::TDomain
    comp::TComp
    params::OrderedDict{Symbol, Parameter}
    cdata
    counter::Int
    lastvalues::Vector{Float64}
    eval::Vector{Float64}
    ipar::Vector{Int}  # handled by Model

    function CompEval(domain::AbstractDomain, comp::AbstractComponent)
        params = getparams(comp)
        (cdata, len) = ceval_data(domain, comp)
        return new{typeof(domain), typeof(comp)}(
            domain, comp, params, cdata, 0,
            fill(NaN, length(params)),
            fill(NaN, len), Vector{Int}())
    end
end


# This is called to `update` to distinguish it from component's `evaluate`.
update(c::CompEval) = update(c, [par.val for par in values(c.params)])
function update(c::CompEval, pvalues::Vector{Float64})
    @assert length(c.params) == length(pvalues)

    # Do we actually need a new evaluation?
    if (any(c.lastvalues .!= pvalues)  ||  (c.counter == 0))
        c.lastvalues .= pvalues
        c.counter += 1
        @assert all(.!isnan.(pvalues))
        evaluate(c, pvalues...)
    end
    return c.eval
end


# ====================================================================
# Component fall back methods
ceval_data(domain::AbstractDomain, comp::AbstractComponent) =
    error("Component " * string(typeof(comp)) * " must implement its own method for `ceval_data`.")

evaluate(c::CompEval{TDomain, TComp}, args...) where {TDomain, TComp} =
    error("Component " * string(TComp) * " must implement its own method for `evaluate`.")


# ====================================================================
# Built-in components
#
include("components/ScalarParam.jl")
include("components/FuncWrap.jl")
include("components/OffsetSlope.jl")
include("components/Gaussian.jl")


# ====================================================================
# Parse a user defined structure or dictionary to extract all
# components
function extract_components(things...; prefix="")
    out = OrderedDict{Symbol, AbstractComponent}()
    for thing in things
        #println()
        #println("Thing $(typeof(thing))  (prefix = $(prefix))")
        if isa(thing, AbstractComponent)
            #println("Adding...")
            out[Symbol(prefix)] = thing
        else
            (length(prefix) > 0)  &&  (prefix *= "_")
            if isa(thing, AbstractDict)
                for (name, v) in thing
                    #println("Dict: Walk through $name :: $(typeof(v))")
                    merge!(out, extract_components(v; prefix=prefix * string(name)))
                end
            elseif isa(thing, Pair)
                #println("Pair: $(thing[1]), $(typeof(thing[2]))")
                if isa(thing[1], Symbol)
                    name = thing[1]
                    v = thing[2]
                    isa(v, Number)  &&  (v = ScalarParam(v))
                    if isa(v, AbstractComponent)
                        merge!(out, extract_components(v; prefix=prefix * string(name)))
                    end
                end
            elseif isstructtype(typeof(thing))
                for name in fieldnames(typeof(thing))
                    v = getfield(thing, name)
                    #println("Structure: Walk through $name :: $(typeof(v))")
                    merge!(out, extract_components(v; prefix=prefix * string(name)))
                end
            end
        end
    end
    return out
end


# ====================================================================
# A model prediction suitable to be compared to experimental data
mutable struct Prediction
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    eval::Vector{Float64}
    reduce_with_dict::Bool
    reducer::Union{Nothing, Function}
    counter::Int

    function Prediction(domain::AbstractDomain, things...;
                        prefix="", reduce=reduce)
        comps = extract_components(things...; prefix=prefix)
        cevals = OrderedDict{Symbol, CompEval}()
        for (name, comp) in comps
            cevals[name] = CompEval(domain, comp)
        end
        out = new(domain, cevals, Vector{Float64}(), false, reduce, 0)
        evaluate(out)  # TODO: is this correct?
        out
    end
end

# Default reducer: add all components
reduce(domain::AbstractDomain, args...) = .+(args...)

# Reduce prediction by combining individual components
function reduce(pred::Prediction)
    if pred.reduce_with_dict
        d = Dict([(cname, ceval.eval) for (cname, ceval) in pred.cevals])
        expr = pred.reducer(pred.domain, d)
    else
        d = [ceval.eval for (cname, ceval) in pred.cevals]
        expr = pred.reducer(pred.domain, d...)
    end
    if length(pred.eval) == 0
        append!(pred.eval, expr)
    else
        pred.eval .= expr
    end
    pred.counter += 1
end

function evaluate(pred::Prediction)
    for (name, ceval) in pred.cevals
        update(ceval)
    end
    reduce(pred)
end


# ====================================================================
# Global model, actually a collection of `Prediction`s.
mutable struct Model
    preds::Vector{Prediction}
    comps::OrderedDict{Symbol, AbstractComponent}
    cfree::OrderedDict{Symbol, Bool}
    params::OrderedDict{NTuple{2, Symbol}, Parameter}
    pvalues::Vector{Float64}
    actual::Vector{Float64}
    buffer::Vector{Float64}
    partransform::Function
end

function Model(v::Vector{Prediction})
    model = Model(v, OrderedDict{Symbol, AbstractComponent}(),
                  OrderedDict{Symbol, Bool}(),
                  OrderedDict{NTuple{2, Symbol}, Parameter}(),
                  Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), default_partransform)
    evaluate(model)
    return model
end

Model(p::Prediction) = Model([p])
Model(args...; kw...) = Model(Prediction(args...; kw...))

function evaluate(model::Model)
    @assert length(model.preds) >= 1

    # Save list of previously free components
    cfree = deepcopy(model.cfree)

    # Collect components and parameters
    empty!(model.comps)
    empty!(model.cfree)
    empty!(model.params)
    for pred in model.preds
        for (cname, ceval) in pred.cevals
            model.comps[cname] = ceval.comp
            model.cfree[cname] = get(cfree, cname, true)
            for (pname, par) in ceval.params
                cpname = (cname, pname)
                model.params[cpname] = par
            end
        end
    end

    # Populate CompEval.ipar and evaluate all predictions
    ndata = 0
    cpnames = keys(model.params)
    for pred in model.preds
        for (cname, ceval) in pred.cevals
            empty!(ceval.ipar)
            for (pname, par) in ceval.params
                cpname = (cname, pname)
                push!(ceval.ipar, findfirst(cpnames .== Ref(cpname)))
            end
            update(ceval)
        end
        reduce(pred)
        ndata += length(pred.eval)
    end

    model.pvalues = [par.val for par in values(model.params)]
    model.actual = deepcopy(model.pvalues)
    model.buffer = Vector{Float64}(undef, ndata)
    quick_evaluate(model)
    return model
end

default_partransform(model::Model, pvalues::Vector{Float64}, actual::Vector{Float64}) = nothing

# This is supposed to be called from `fit!`, not by user
function quick_evaluate(model::Model)
    model.actual .= model.pvalues  # copy all values by default
    model.partransform(model, model.pvalues, model.actual)

    for pred in model.preds
        for (cname, ceval) in pred.cevals
            update(ceval, model.actual[ceval.ipar])
        end
    end
    for pred in model.preds
        reduce(pred)
    end
    nothing
end


function Base.push!(m::Model, p::Prediction)
    push!(m.preds, p)
    evaluate(model)
    return model
end

Base.getindex(m::Model, i::Int) = m.preds[i].eval
Base.getindex(m::Model, cname::Symbol) = m.comps[cname]
Base.getindex(m::Model, cname::Symbol, pname::Symbol) = m.params[(cname, pname)]

parindex(model::Model, cname::Symbol, pname::Symbol) =
    findfirst(keys(model.params) .== Ref((cname, pname)))

function freeze(model::Model, cname::Symbol)
    @assert cname in keys(model.cfree) "Component $c is not defined"
    model.cfree[cname] = false
    model
end

function thaw(model::Model, cname::Symbol)
    @assert cname in keys(model.cfree) "Component $c is not defined"
    model.cfree[cname] = true
    model
end

# ====================================================================
# Fit results
#
struct BestFitPar
    val::Float64
    unc::Float64
    free::Bool
    actual::Float64  # value after transformation
end

struct BestFitResult
    comps::OrderedDict{Symbol, OrderedDict{Symbol, Union{BestFitPar, Vector{BestFitPar}}}}
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:Optimal, :NonOptimal, :Warn, :Error
    log10testprob::Float64
    elapsed::Float64
end

Base.getindex(res::BestFitResult, cname::Symbol) = res.comps[cname]

# ====================================================================
function data1D(model::Model, data::Vector{T}) where T<:AbstractMeasures
    out = Vector{Measures_1D}()
    for i in 1:length(model.preds)
        pred = model.preds[i]
        @assert(length(data[i]) == length(pred.eval),
                "Length of dataset $i do not match corresponding model prediction.")
        push!(out, flatten(data[i], pred.domain))
    end
    return out
end


function residuals1d(model::Model, data1d::Vector{Measures_1D})
    c1 = 1
    for i in 1:length(model.preds)
        pred = model.preds[i]
        c2 = c1 + length(pred.eval) - 1
        model.buffer[c1:c2] .= ((pred.eval .- data1d[i].val) ./ data1d[i].unc)
        c1 = c2 + 1
    end
    return model.buffer
end



# ====================================================================
abstract type AbstractMinimizer end

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


fit!(model::Model, data::T; kw...) where T<:AbstractMeasures =
    fit!(model, [data]; kw...)

function fit!(model::Model, data::Vector{T};
              minimizer=lsqfit()) where T<:AbstractMeasures
    elapsedTime = Base.time_ns()
    evaluate(model)

    free = Vector{Bool}()
    for (cpname, par) in model.params
        push!(free, par.free  &&  model.cfree[cpname[1]])
    end
    ifree = findall(free)
    @assert length(ifree) > 0 "No free parameter in the model"

    # Flatten empirical data
    data1d = data1D(model, data)

    # Evaluate normalized residuals starting from free parameter values
    function pval2resid(pvalues_free::Vector{Float64})
        model.pvalues[ifree] .= pvalues_free  # update parameter values
        quick_evaluate(model)
        return residuals1d(model, data1d)
    end

    (status, best_val, best_unc) = minimize(minimizer, pval2resid,
                                            collect(values(model.params))[ifree])

    model.pvalues[ifree] .= best_val
    setfield!.(values(model.params), :val, model.pvalues)
    uncerts = fill(NaN, length(model.pvalues))
    uncerts[ifree] .= best_unc

    # Prepare output
    quick_evaluate(model)  # ensure best fit values are used
    comp = OrderedDict{Symbol, OrderedDict{Symbol, BestFitPar}}()
    for cname in keys(model.comps)
        comp[cname] = OrderedDict{Symbol, BestFitPar}()
    end
    i = 1
    for (cpname, par) in model.params
        cname = cpname[1]
        pname = cpname[2]
        comp[cname][pname] = BestFitPar(model.pvalues[i], uncerts[i],
                                        (i in ifree), model.actual[i])
        i += 1
    end
    cost = sum(abs2, model.buffer)
    dof = length(model.buffer) - length(ifree)
    result = BestFitResult(comp, length(model.buffer), dof, cost, status,
                           logccdf(Chisq(dof), cost) * log10(exp(1)),
                           float(Base.time_ns() - elapsedTime) / 1.e9)
    return result
end

include("show.jl")

end
