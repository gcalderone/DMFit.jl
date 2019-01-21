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


# --------------------------------------------------------------------
# Show methods

printmain(io::IO, args...)                  =                    printstyled(io, args..., "\n"; bold=true, color=:light_blue)
printerr( io::IO, args...)                  =                    printstyled(io, args..., "\n"; color=:red, bold=true)
function printtype(io::IO, args...)         ; print(io, "    "); printstyled(io, args..., "\n"; bold=true); end
function printhead(io::IO, args...)         ; print(io, "    "); printstyled(io, args..., "\n"; color=:underline); end
function printcell(io::IO, args...; u=false); print(io, "    "); printstyled(io, args..., "\n"; color=(u  ?  :underline  :  :default)); end


function show(io::IO, dom::AbstractCartesianDomain)
    printtype(io, typeof(dom), "  length: ", length(dom), "")
    s = @sprintf("%6s │ %8s │ %10s │ %10s │ %10s │ %10s",
                 "Dim.", "Size", "Min val", "Max val", "Min step", "Max step")
    printhead(io, s)

    for i in 1:ndims(dom)
        a = dom[i]
        b = 0
        if length(a) > 1
            b = a .- circshift(a, 1)
            b = b[2:end]
        end
        s = @sprintf("%6d │ %8d │ %10.4g │ %10.4g │ %10.4g │ %10.4g",
                     i, length(a),
                     minimum(a), maximum(a),
                     minimum(b), maximum(b))
        printcell(io, s)
    end
end


function show(io::IO, dom::AbstractLinearDomain)
    printtype(io, typeof(dom), " length: ", length(dom), "")
    s = @sprintf("%6s │ %10s │ %10s",
                 "Dim.", "Min val", "Max val")
    printhead(io, s)
    
    for i in 1:ndims(dom)
        s = @sprintf("%6d │ %10.4g │ %10.4g",
                     i, getaxismin(dom, i), getaxismax(dom, i))
        printcell(io, s)
    end
end


# Special case for Domain_1D: treat it as a Cartesian domain, despite it is a Linear one.
function show(io::IO, dom::Domain_1D)
    printtype(io, typeof(dom), " length: ", length(dom), "")
    s = @sprintf("%6s │ %8s │ %10s │ %10s │ %10s │ %10s",
                 "Dim.", "Size", "Min val", "Max val", "Min step", "Max step")
    printhead(io, s)

    a = dom[1]
    b = 0
    if length(a) > 1
        b = a .- circshift(a, 1)
        b = b[2:end]
    end
    s = @sprintf("%6d │ %8d │ %10.4g │ %10.4g │ %10.4g │ %10.4g",
                 1, length(a),
                 minimum(a), maximum(a),
                 minimum(b), maximum(b))
    printcell(io, s)
end


# --------------------------------------------------------------------
function show(io::IO, data::AbstractData)
    printtype(io, typeof(data), "   length: ", (length(data.measure)))
    printhead(io, @sprintf("%8s │ %10s │ %10s │ %10s │ %10s │ %10s",
                           "", "Min", "Max", "Mean", "Median", "Stddev."))

    nonFinite = Vector{String}()
    names = fieldnames(typeof(data))
    for name in names
        a = getfield(data, name)

        nan = length(findall(isnan.(a)))
        inf = length(findall(isinf.(a)))

        if nan > 0  || inf > 0
            push!(nonFinite, @sprintf("%8s │ NaN: %-10d   Inf: %-10d",
                                      string(name), nan, inf))
            a = a[findall(isfinite.(a))]
        end

        s = @sprintf("%8s │ %10.4g │ %10.4g │ %10.4g │ %10.4g │ %10.4g",
                     string(name),
                     minimum(a), maximum(a),
                     mean(a), median(a), std(a))
        printcell(io, s)
    end
    
    if length(nonFinite) > 0
        println(io)
        for s in nonFinite
            printerr(s)
        end
    end
end




# ====================================================================
# Parameter structure and show method for Component objects
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
    

function show(io::IO, comp::AbstractComponent; header=true, count=0, cname="")
    (header)  &&  (printtype(io, typeof(comp)))
    if count == 0
        printhead(io, @sprintf "%5s │ %20s │ %10s │ %10s │ %10s │ %10s │ %s"  "#" "Component" "Param." "Value" "Low" "High" "Notes")
    end

    localcount = 0; lastcount = length(getparams(comp))
    for (pname, wparam) in getparams(comp)
        localcount += 1
        par = wparam.par
        note = ""
        (par.fixed)  &&  (note *= "FIXED")
        (par.expr != "")  &&  (note *= " expr=" * par.expr)
        count += 1
        s = @sprintf("%5d │ %20s │ %10s │ %10.3g │ %10.3g │ %10.3g │ %s",
                     count, cname,
                     (wparam.index >= 1  ?  Symbol(wparam.pname, "[", wparam.index, "]")  :  pname),
                     par.val, par.low, par.high, note)
        
        printcell(io, s, u=(localcount == lastcount))
    end
    return count
end


# ====================================================================
# CompEvaluation, CompiledExpression and Model structure, and
# associated show method
#
mutable struct CompEvaluation
    counter::Int
    cdata::AbstractComponentData
    lastParams::Vector{Float64}
    result::Vector{Float64}
end


mutable struct CompiledExpression
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
    compiled::Vector{CompiledExpression}
    Model() = new(OrderedDict{Symbol, AbstractComponent}(),
                  Vector{CompiledExpression}())
end


show(io::IO, mime::MIME"text/plain", model::Model) = show(io, model)
function show(io::IO, model::Model)
    printmain(io, @sprintf "Components:")
    compcount(model) != 0  || (return nothing)

    printhead(io, @sprintf "%5s │ %20s │ %s"  "#" "Component" "Type")

    count = 0
    for (cname, comp) in components(model)
        count += 1
        ctype = split(string(typeof(comp)), ".")
        (ctype[1] == "DataFitting")  &&   (ctype = ctype[2:end])
        ctype = join(ctype, ".")

        s = @sprintf "%5d │ %20s │ %s" count string(cname) ctype
        printcell(io, s)
    end
    println(io)

    printmain(io, "Parameters:")
    count = 0
    for (cname, comp) in components(model)
        count = show(io, comp, cname=string(cname), count=count, header=false)
    end
    println(io)
    
    if length(compiled(model)) == 0
        printmain(io, "Total expressions: 0")
        return nothing
    end
    
    printmain(io, "Domain(s):")
    for ii in 1:length(compiled(model))
        ce = compiled(model, ii)
        printmain(io, @sprintf "#%d:  " ii)
        show(io, ce.domain)
        println(io)
    end

    countexpr = 0
    printmain(io, "Expression(s): ")
    printhead(io, @sprintf "%3s │ %10s │ %7s │ %10s │ %10s │ %10s │ %10s │ %10s │ " "#" "Component" "Counter" "Min" "Max" "Mean" "NaN" "Inf")

    for ii in 1:length(compiled(model))
        ce = compiled(model, ii)

        for jj in 1:length(ce.compevals)
            cname = ce.compnames[jj]
            ceval = ce.compevals[jj]

            result = ceval.result
            v = view(result, findall(isfinite.(result)))
            nan = length(findall(isnan.(result)))
            inf = length(findall(isinf.(result)))
            printcell(io, @sprintf("%3d │ %10s │ %7d │ %10.3g │ %10.3g │ %10.3g │ %10d │ %10d │ ",
                                   ii, cname, ceval.counter,
                                   minimum(v), maximum(v), mean(v), nan, inf))
        end

        localcount = 0; lastcount = length(ce.exprs)
        for jj in 1:length(ce.exprs)
            localcount += 1
            result = ce.results[jj]
            v = view(result, findall(isfinite.(result)))
            nan = length(findall(isnan.(result)))
            inf = length(findall(isinf.(result)))
            printcell(io, u=(localcount == lastcount),
                      @sprintf("%3d │ %10s │ %7d │ %10.3g │ %10.3g │ %10.3g │ %10d │ %10d │ %s",
                               ii, "Expr #"*string(jj), ce.counter,
                               minimum(v), maximum(v), mean(v), nan, inf, ce.exprs[jj]))
            countexpr += 1
        end
    end
    println(io)
    printmain(io, "Total expressions: ", countexpr)
end


# ====================================================================
# FitParameter and FitResult structures, and associated show method.
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


function show(io::IO, comp::BestFitComp; count=0, cname="")
    if count == 0
        printhead(io, @sprintf "%5s │ %20s │ %10s │ %10s │ %10s │ %10s"  "#" "Component" "Param." "Value" "Uncert." "Rel. unc. (%)")
    end
    localcount = 0;  lastcount = length(getfield(comp, :params))
    for (pname, params) in getfield(comp, :params)
        localcount += 1
        if typeof(params) == Vector{BestFitParam}
            for ii in 1:length(params)
                count += 1
                par = params[ii]
                spname = string(pname) * "[" * string(ii) * "]"
                printcell(io, u=((localcount == lastcount)  &&  (ii == length(params))),
                                 @sprintf("%5d │ %20s │ %10s │ %10.4g │ %10.4g │ %10.2g", count, cname,
                                          spname, par.val, par.unc, par.unc/par.val*100.))
            end
        else
            count += 1
            par = params
            spname = string(pname)
            s = @sprintf("%5d │ %20s │ %10s │ %10.4g │ %10.4g │ %10.2g", count, cname,
                         spname, par.val, par.unc, par.unc/par.val*100.)
            printcell(io, s, u=(localcount == lastcount))
        end
    end
    return count
end


function show(io::IO, f::FitResult)
    printmain(io, @sprintf "Best fit values:")

    count = 0
    for (cname, comp) in getfield(f.bestfit, :comp)
        count = show(io, comp, count=count, cname=string(cname))
    end

    println(io)
    printmain(io, "Summary:")
    printcell(io, @sprintf("#Data  : %10d              Cost: %10.5g", f.ndata, f.cost))
    printcell(io, @sprintf("#Param : %10d              DOF : %10d", f.ndata-f.dof, f.dof))
    printcell(io, @sprintf("Elapsed: %10.4g s            Red.: %10.4g", f.elapsed, f.cost / f.dof))
    printstyled(io, "    Status :  ", bold=true)
    if f.status == :Optimal
        printstyled(color=:green, io, "Optimal", bold=true)
    elseif f.status == :NonOptimal
        printstyled(color=:yellow, io, "non-Optimal, see fitter output", bold=true)
    elseif f.status == :Warn
        printstyled(color=:yellow, io, "Warning, see fitter output", bold=true)
    elseif f.status == :Error
        printstyled(color=:red, io, "Error, see fitter output", bold=true)
    else
        printstyled(color=:magenta, io, "Unknown (" * string(f.status) * "), see fitter output", bold=true)
    end
    println(io)
end



# ====================================================================
# Built-in components: SimpleParam and FuncWrapper
#

# --------------------------------------------------------------------
# SimpleParam
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



# --------------------------------------------------------------------
# FuncWrap
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
