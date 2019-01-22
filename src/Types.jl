# ====================================================================
# Abstract types
#
abstract type AbstractDomain end
abstract type AbstractLinearDomain    <: AbstractDomain end
abstract type AbstractCartesianDomain <: AbstractDomain end

abstract type AbstractData end
abstract type AbstractMeasures <: AbstractData end
abstract type AbstractCounts   <: AbstractData end

abstract type AbstractComponent end
abstract type AbstractComponentData end

abstract type AbstractMinimizer end


# ====================================================================
# Define domain, data and associated methods for ndim=1, 2 and 3
#

# The following is used in `code_ndim` macro, hence it must be declared here
mutable struct FuncWrap_cdata <: AbstractComponentData
    func::Function
end


macro code_ndim(ndim::Int)
    @assert ndim >= 1 "Number of dimensions must be >= 1"
    out = Expr(:block)

    # Structure def.
    if ndim > 1
        name = Symbol(:CartesianDomain_, ndim, :D)
        push!(out.args, :(
            struct $name <: AbstractCartesianDomain
                axis::Vector{Vector{Float64}}
                size::NTuple{$(ndim), Int}
                index::Vector{Int}
            end;
            Base.ndims(dom::$name) = $ndim;
            Base.size(dom::$name) = dom.size;
            Base.size(dom::$name, dim::Int) = dom.size[dim];
            Base.getindex(dom::$name, dim::Int) = dom.axis[dim];
        ))

        # Constructors
        a = ["d$(i)::AbstractVector{Float64}" for i in 1:ndim]
        b = ["deepcopy(d$(i))" for i in 1:ndim]
        c = ["length(d$(i))" for i in 1:ndim]
        s = "CartesianDomain(" * join(a, ", ") * "; index=1:" * join(c, " * ") *
            ") = CartesianDomain_$(ndim)D([" * join(b, ", ") * "], (" * join(c, ", ") * ",), index)"
        push!(out.args, Meta.parse(s))

        a = ["n$(i)::Int" for i in 1:ndim]
        b = ["one(Float64):one(Float64):n$(i)" for i in 1:ndim]
        c = ["n$(i)" for i in 1:ndim]
        s = "CartesianDomain(" * join(a, ", ") * "; index=1:" * join(c, " * ") *
            ") = CartesianDomain(" * join(b, ", ") * ", index=index)"
        push!(out.args, Meta.parse(s))
    end

    # Structure def.
    name = Symbol(:Domain_, ndim, :D)
    tmp = (ndim > 1  ?  2  :  1)
    push!(out.args, :(
        struct $name <: AbstractLinearDomain
            axis::Array{Float64, $tmp}
            vmin::Vector{Float64}
            vmax::Vector{Float64}
            length::Int
        end;
        Base.ndims(dom::$name) = $ndim;
    ))

    # getindex
    if ndim == 1
        push!(out.args, :(Base.getindex(dom::$name, dim::Int) = dom.axis))
    else
        push!(out.args, :(Base.getindex(dom::$name, dim::Int) = dom.axis[dim,:]))
    end

    # Constructors
    a = ["d$(i)::AbstractVector{Float64}" for i in 1:ndim]
    b = ["deepcopy(d$(i))" for i in 1:ndim]
    c = ["d$(i)" for i in 1:ndim]
    s = "function Domain(" * join(a, ", ") * ")\n"
    if ndim > 1
        s *= "  @assert " * join("length(" .* c .* ")", " == ") * " \"Arrays must have same length\" \n"
    end
    s *= "  return Domain_$(ndim)D(" *
        (ndim == 1  ?  "d1"  :  "[" * join(b, " ") * "]'") * ", " *
        "[" * join("minimum(" .* c .* ")", ", ") * "], " *
        "[" * join("maximum(" .* c .* ")", ", ") * "], " *
        "length(d1))\n"
    s *= "end"
    push!(out.args, Meta.parse(s))

    a = ["n$(i)::Int" for i in 1:ndim]
    b = ["one(Float64):one(Float64):n$(i)" for i in 1:ndim]
    s = "Domain(" * join(a, ", ") * ") = Domain(" * join(b, ", ") * ")"
    push!(out.args, Meta.parse(s))

    # flatten
    if ndim > 1
        s = "function flatten(dom::CartesianDomain_$(ndim)D)::Domain_$(ndim)D\n"
        s *= "  out = Matrix{Float64}(undef, $ndim, length(dom))\n"
        s *= "  for i in 1:length(dom)\n"
        a = ["d$(i)" for i in 1:ndim]
        s *= "  (" * join(a, ", ") * ") = Tuple(CartesianIndices(size(dom))[dom.index[i]])\n"
        for i in 1:ndim
            s *= "  out[$i, i] = dom.axis[$i][d$(i)]\n"
        end
        s *= "  end\n"
        a = ["out[$(i), :]" for i in 1:ndim]
        s *= "  return Domain(" * join(a, ", ") * ")\n"
        s *= "end\n"
        push!(out.args, Meta.parse(s))
    end

    name = Symbol(:Measures_, ndim, :D)
    push!(out.args, :(
        struct $name <: AbstractMeasures
            measure::Array{Float64, $(ndim)}
            uncert::Array{Float64, $(ndim)}
        end;
        Base.ndims(dom::$name) = $ndim;
    ))

    push!(out.args, :(
        function Measures(measure::Array{Float64, $ndim}, uncert::Array{Float64, $(ndim)})
            @assert length(measure) == length(uncert) "Measure and uncertainty arrays must have same size"
            $(name)(measure, uncert)
        end;
    ))

    push!(out.args, :(
        function Measures(measure::Array{Float64, $ndim}, uncert::Float64=one(Float64))
            u = similar(measure)
            fill!(u, uncert)
            $(name)(measure, u)
        end;
    ))

    name = Symbol(:Counts_, ndim, :D)
    push!(out.args, :(
        struct $name <: AbstractCounts
            measure::Array{Int, $(ndim)}
        end;
        Base.ndims(dom::$name) = $ndim;
    ))

    push!(out.args, :(
        function Counts(measure::Array{Int, $ndim})
            $(name)(measure)
        end
    ))

    s = Vector{String}()
    push!(s, "function evaluate!(output::Vector{Float64}, domain::Domain_$(ndim)D, compdata::FuncWrap_cdata, params...)")
    push!(s, "  output .= compdata.func(" * join("domain[" .* string.(collect(1:ndim)) .* "]", ", ") * ", params...)")
    push!(s, "end")
    push!(out.args, Meta.parse(join(s, "\n")))

    return esc(out)
end

@code_ndim 1
@code_ndim 2

# The following methods do not require a macro to be implemented
getaxismin(dom::AbstractLinearDomain, dim::Int) = dom.vmin[dim]
getaxismax(dom::AbstractLinearDomain, dim::Int) = dom.vmax[dim]
getaxisextrema(dom::AbstractLinearDomain, dim::Int) = (dom.vmin[dim], dom.vmax[dim])
length(dom::AbstractLinearDomain) = dom.length
length(dom::AbstractCartesianDomain) = length(dom.index)
length(data::AbstractData) = length(data.measure)
size(data::AbstractData) = size(data.measure)
size(data::AbstractData, dim::Int) = size(data.measure)[dim]

# Methods to "flatten" a multidimensional object <: AbstractData into a 1D one
flatten(data::AbstractMeasures, dom::AbstractCartesianDomain)::Measures_1D = Measures_1D(data.measure[dom.index], data.uncert[dom.index])
flatten(data::AbstractCounts, dom::AbstractCartesianDomain)::Counts_1D = Counts_1D(data.measure[dom.index])

"""
# reshape

Reshape an array according to the size of a CartesianDomain object
"""
function reshape(array::AbstractArray, dom::AbstractCartesianDomain)
    @assert length(array) == length(dom) "Domain and array must have the same length"
    out = Array{Float64}(undef, size(dom)...)
    out .= NaN
    out[dom.index] .= array
    return out
end


# ====================================================================
# Parameter and WrapParameter structure
#
mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    step::Float64
    fixed::Bool               # true = fixed; false = free
    expr::String
    Parameter(value::Number) = new(float(value), -Inf, +Inf, NaN, false, "")
end


mutable struct WrapParameter
    pname::Symbol
    index::Int
    par::Parameter
end


# ====================================================================
# CompEvaluation, Instrument and Model structure
#
mutable struct CompEvaluation
    counter::Int
    cdata::AbstractComponentData
    lastParams::Vector{Float64}
    result::Vector{Float64}
end


mutable struct Instrument
    label::String
    code::String
    funct::Function
    counter::Int

    domain::AbstractDomain
    exprs::Vector{Expr}
    results::Vector{Vector{Float64}}

    compnames::Vector{Symbol}  # only involved components
    compevals::Vector{CompEvaluation}
end


mutable struct Model <: AbstractDict{Symbol, AbstractComponent}
    comp::OrderedDict{Symbol, AbstractComponent}
    instruments::Vector{Instrument}
    Model() = new(OrderedDict{Symbol, AbstractComponent}(),
                  Vector{Instrument}())
end


# ====================================================================
# Fit results
#
struct BestFitParam
    val::Float64
    unc::Float64
end

struct BestFitComp
    params::OrderedDict{Symbol, Union{BestFitParam, Vector{BestFitParam}}}
end

struct BestFit
    comp::OrderedDict{Symbol, BestFitComp}
end

struct FitResult
    fitter::AbstractMinimizer
    bestfit::BestFit
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:Optimal, :NonOptimal, :Warn, :Error
    elapsed::Float64
end


# ____________________________________________________________________
# Built-in component: SimpleParam
mutable struct SimpleParam_cdata <: AbstractComponentData; end

mutable struct SimpleParam <: AbstractComponent
    val::Parameter
    SimpleParam(val::Number) = new(Parameter(val))
end

compdata(domain::AbstractDomain, comp::SimpleParam) = SimpleParam_cdata()
function evaluate!(output::AbstractArray{Float64}, domain::AbstractDomain,
                   compdata::SimpleParam_cdata, val)
    output .= val
    return output
end


# ____________________________________________________________________
# Built-in component: FuncWrap
mutable struct FuncWrap <: AbstractComponent
    func::Function
    p::Vector{Parameter}
end

function FuncWrap(func::Function, args...)
    params= Vector{Parameter}()
    for i in 1:length(args)
        push!(params, Parameter(args[i]))
    end
    return FuncWrap(func, params)
end

compdata(domain::DataFitting.AbstractDomain, comp::FuncWrap) = FuncWrap_cdata(comp.func)
