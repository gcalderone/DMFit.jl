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


macro code_ndim(ndim::Int, importf=true)
    prefix = (importf  ?  "ModelFit."  :  "")
    
    @assert ndim >= 1 "Number of dimensions must be >= 1"
    out = Expr(:block)

    if importf
        s = Vector{String}()
        push!(s, "import ModelFit.Domain, ModelFit.flatten, ModelFit.Measures, ModelFit.Counts, ModelFit.evaluate!")
        push!(out.args, Meta_parse(join(s, "\n")))
    end

    # Structure def.
    if ndim > 1
        s = Vector{String}()
        if importf
            push!(s, "import ModelFit.CartesianDomain")
            push!(out.args, Meta_parse(join(s, "\n")))
        end
        
        s = Vector{String}()
        push!(s, "struct CartesianDomain_$(ndim)D <: $(prefix)AbstractCartesianDomain")
        push!(s, "  d::Vector{Vector{Float64}}")
        push!(s, "  size::NTuple{$(ndim), Int}")
        push!(s, "  index::Vector{Int}")
        push!(s, "end")
        push!(out.args, Meta_parse(join(s, "\n")))

        # size
        s = "Base.size(dom::CartesianDomain_$(ndim)D) = dom.size"
        push!(out.args, Meta_parse(s))
        s = "Base.size(dom::CartesianDomain_$(ndim)D, dim::Int) = dom.size[dim]"
        push!(out.args, Meta_parse(s))

        # ndims
        s = "Base.ndims(dom::CartesianDomain_$(ndim)D) = $ndim"
        push!(out.args, Meta_parse(s))

        # getindex
        s = "Base.getindex(dom::CartesianDomain_$(ndim)D, dim::Int) = dom.d[dim]"
        push!(out.args, Meta_parse(s))

        # Constructors
        a = Vector{String}()
        for i in 1:ndim
            push!(a, "d$(i)::AbstractVector{Float64}")
        end
        b = Vector{String}()
        for i in 1:ndim
            push!(b, "deepcopy(d$(i))")
        end
        c = Vector{String}()
        for i in 1:ndim
            push!(c, "length(d$(i))")
        end
        s = "CartesianDomain(" * join(a, ", ") * "; index=1:" * join(c, " * ") * ") = CartesianDomain_$(ndim)D([" * join(b, ", ") * "], (" * join(c, ", ") * ",), index)"
        push!(out.args, Meta_parse(s))


        a = Vector{String}()
        for i in 1:ndim
            push!(a, "n$(i)::Int")
        end
        b = Vector{String}()
        for i in 1:ndim
            push!(b, "one(Float64):one(Float64):n$(i)")
        end
        c = Vector{String}()
        for i in 1:ndim
            push!(c, "n$(i)")
        end
        s = "CartesianDomain(" * join(a, ", ") * "; index=1:" * join(c, " * ") * ") = CartesianDomain(" * join(b, ", ") * ", index=index)"
        push!(out.args, Meta_parse(s))
    end

    # Structure def.
    s = Vector{String}()
    push!(s, "struct Domain_$(ndim)D <: $(prefix)AbstractLinearDomain")
    if ndim == 1
        push!(s, "  d::Vector{Float64}")
    else
        push!(s, "  d::Matrix{Float64}")
    end
    push!(s, "  vmin::Vector{Float64}")
    push!(s, "  vmax::Vector{Float64}")
    push!(s, "  length::Int")
    push!(s, "end")
    push!(out.args, Meta_parse(join(s, "\n")))

    # ndims
    s = "Base.ndims(dom::Domain_$(ndim)D) = $ndim"
    push!(out.args, Meta_parse(s))

    # getindex
    if ndim == 1
        s = "Base.getindex(dom::Domain_1D, dim::Int) = dom.d"
    else
        s = "Base.getindex(dom::Domain_$(ndim)D, dim::Int) = dom.d[dim,:]"
    end
    push!(out.args, Meta_parse(s))

    # Constructors
    a = Vector{String}()
    for i in 1:ndim
        push!(a, "d$(i)::AbstractVector{Float64}")
    end
    b = Vector{String}()
    for i in 1:ndim
        push!(b, "deepcopy(d$(i))")
    end
    c = Vector{String}()
    for i in 1:ndim
        push!(c, "d$(i)")
    end
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
    push!(out.args, Meta_parse(s))

    a = Vector{String}()
    for i in 1:ndim
        push!(a, "n$(i)::Int")
    end
    b = Vector{String}()
    for i in 1:ndim
        push!(b, "one(Float64):one(Float64):n$(i)")
    end
    s = "Domain(" * join(a, ", ") * ") = Domain(" * join(b, ", ") * ")"
    push!(out.args, Meta_parse(s))


    # flatten
    if ndim > 1
        s = "function flatten(dom::CartesianDomain_$(ndim)D)::Domain_$(ndim)D\n"
        if oldver()
            s *= "  out = Matrix{Float64}(       $ndim, length(dom))\n"
        else
            s *= "  out = Matrix{Float64}(undef, $ndim, length(dom))\n"
        end
        s *= "  for i in 1:length(dom)\n"
        a = Vector{String}()
        for i in 1:ndim
            push!(a, "d$(i)")
        end
        if oldver()
            s *= "  (" * join(a, ", ") * ") = ind2sub(size(dom), dom.index[i])\n"
        else
            s *= "  (" * join(a, ", ") * ") = Tuple(CartesianIndices(size(dom))[dom.index[i]])\n"
        end
        for i in 1:ndim
            s *= "  out[$i, i] = dom.d[$i][d$(i)]\n"
        end
        s *= "  end\n"
        a = Vector{String}()
        for i in 1:ndim
            push!(a, "out[$(i), :]")
        end
        s *= "  return Domain(" * join(a, ", ") * ")\n"

        s *= "end\n"
        push!(out.args, Meta_parse(s))
    end

    s = Vector{String}()
    push!(s, "struct Measures_$(ndim)D <: $(prefix)AbstractMeasures")
    push!(s, "  measure::Array{Float64, $(ndim)}")
    push!(s, "  uncert::Array{Float64, $(ndim)}")
    push!(s, "end")
    push!(out.args, Meta_parse(join(s, "\n")))

    s = Vector{String}()
    push!(s, "  function Measures(measure::Array{Float64, $(ndim)}, uncert::Array{Float64, $(ndim)})")
    push!(s, "      @assert length(measure) == length(uncert) \"Measure and uncertainty arrays must have same size\"")
    push!(s, "      $(prefix)Measures_$(ndim)D(measure, uncert)")
    push!(s, "  end")
    push!(out.args, Meta_parse(join(s, "\n")))

    s = Vector{String}()
    push!(s, "function Measures(measure::Array{Float64, $(ndim)}, uncert::Float64=one(Float64))")
    push!(s, "    e = similar(measure)")
    push!(s, "    fill!(e, uncert)")
    push!(s, "    Measures_$(ndim)D(measure, e)")
    push!(s, "end")
    push!(out.args, Meta_parse(join(s, "\n")))

    s = Vector{String}()
    push!(s, "struct Counts_$(ndim)D <: $(prefix)AbstractCounts")
    push!(s, "  measure::Array{Int, $(ndim)}")
    push!(s, "end")
    push!(out.args, Meta_parse(join(s, "\n")))

    s = Vector{String}()
    push!(s, "function Counts(measure::Array{Int, $(ndim)})")
    push!(s, "    Counts_$(ndim)D(measure)")
    push!(s, "end")
    push!(out.args, Meta_parse(join(s, "\n")))

    s = "Base.ndims(dom::Measures_$(ndim)D) = $ndim"
    push!(out.args, Meta_parse(s))
    s = "Base.ndims(dom::Counts_$(ndim)D) = $ndim"
    push!(out.args, Meta_parse(s))

    s = Vector{String}()
    push!(s, "function evaluate!(output::Vector{Float64}, domain::Domain_$(ndim)D, compdata::$(prefix)FuncWrap_cdata, params...)")
    push!(s, "  output .= compdata.func(" * join("domain[" .* string.(collect(1:ndim)) .* "]", ", ") * ", params...)")
    push!(s, "end")
    push!(out.args, Meta_parse(join(s, "\n")))

    return esc(out)
end

@code_ndim 1 false
@code_ndim 2 false

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

Reshape an array according to the size of a CartedsianDomain object
"""
function reshape(array::AbstractArray, dom::AbstractCartesianDomain)
    @assert length(array) == length(dom) "Domain and array must have the same length"
    if oldver()
        out = Array{Float64}(       size(dom)...)
    else
        out = Array{Float64}(undef, size(dom)...)
    end
    out .= NaN
    out[dom.index] .= array
    return out
end


# --------------------------------------------------------------------
# Show methods
function show(stream::IO, dom::AbstractCartesianDomain)
    printstyled(color=:default, stream, bold=true, typeof(dom), "  length: ", length(dom), "\n")
    s = @sprintf("%6s|%8s|%10s|%10s|%10s|%10s",
                 "Dim.", "Size", "Min val", "Max val", "Min step", "Max step")
    println(stream, s)

    for i in 1:ndims(dom)
        a = dom[i]
        b = 0
        if length(a) > 1
            b = a .- circshift(a, 1)
            b = b[2:end]
        end
        s = @sprintf("%6d|%8d|%10.4g|%10.4g|%10.4g|%10.4g",
                     i, length(a),
                     minimum(a), maximum(a),
                     minimum(b), maximum(b))
        println(stream, s)
    end
end


function show(stream::IO, dom::AbstractLinearDomain)
    printstyled(color=:default, stream, bold=true, typeof(dom), " length: ", length(dom), "\n")
    s = @sprintf("%6s|%10s|%10s",
                 "Dim.", "Min val", "Max val")
    println(stream, s)

    for i in 1:ndims(dom)
        s = @sprintf("%6d|%10.4g|%10.4g",
                     i, getaxismin(dom, i), getaxismax(dom, i))

        println(stream, s)
    end
end


# Special case for Domain_1D: treat it as a Cartesian domain, despite it is a Linear one.
function show(stream::IO, dom::Domain_1D)
    printstyled(color=:default, stream, bold=true, typeof(dom), " length: ", length(dom), "\n")
    s = @sprintf("%6s|%8s|%10s|%10s|%10s|%10s",
                 "Dim.", "Size", "Min val", "Max val", "Min step", "Max step")
    println(stream, s)

    a = dom[1]
    b = 0
    if length(a) > 1
        b = a .- circshift(a, 1)
        b = b[2:end]
    end
    s = @sprintf("%6d|%8d|%10.4g|%10.4g|%10.4g|%10.4g",
                 1, length(a),
                 minimum(a), maximum(a),
                 minimum(b), maximum(b))
    println(stream, s)
end



function show(stream::IO, data::AbstractData)
    printstyled(color=:default, stream, bold=true, typeof(data))
    println(stream, "   length: ", (length(data.measure)))
    s = @sprintf("%8s | %10s | %10s | %10s | %10s | %10s",
                 "", "Min", "Max", "Mean", "Median", "Stddev.")
    println(stream, s)

    nonFinite = Vector{String}()
    names = fieldnames(typeof(data))
    for name in names
        a = getfield(data, name)

        nan = length(findall(isnan.(a)))
        inf = length(findall(isinf.(a)))

        if nan > 0  || inf > 0
            push!(nonFinite, @sprintf("%8s | NaN: %-10d   Inf: %-10d\n",
                                      string(name), nan, inf))

            a = a[findall(isfinite.(a))]
        end

        s = @sprintf("%8s | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g",
                     string(name),
                     minimum(a), maximum(a),
                     mean(a), median(a), std(a))
        println(stream, s)
    end

    if length(nonFinite) > 0
        println(stream)
        for s in nonFinite
            printstyled(color=:red, stream, bold=true, s)
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


function show(stream::IO, comp::AbstractComponent; color=:default, header=true, count=0, cname="")
    (header)  &&  (printstyled(color=:default, stream, bold=true, typeof(comp)); println(stream))
    if count == 0
        s = @sprintf "%5s|%20s|%10s|%10s|%10s|%10s|%s\n"  "#" "Component" "Param." "Value" "Low" "High" "Notes"
        printstyled(color=:default, stream, s)
    end

    extraFields = false
    (wnames, params) = getparams(comp)
    for i in 1:length(params)
        par = params[i]
        wname = wnames[i]

        note = ""
        (par.fixed)  &&  (note *= "FIXED")
        (par.expr != "")  &&  (note *= " expr=" * par.expr)

        count += 1
        s = @sprintf("%5d|%20s|%10s|%10.3g|%10.3g|%10.3g|%s\n",
                     count, cname, wname,
                     par.val, par.low, par.high, note)
        printstyled(color=color, stream, s)
    end

    if extraFields  &&  header
        println(stream)
        println(stream, "Extra fields:")
        for pname in fieldnames(typeof(comp))
            if fieldtype(typeof(comp), pname) == Parameter
                continue
            end
            println(stream, string(pname), " : ", string(fieldtype(typeof(comp), pname)))
        end
    end
    return count
end


# ====================================================================
# ComponentEvaluation, CompiledExpression and Model structure, and
# associated show method
#
mutable struct ComponentEvaluation
    cdata::AbstractComponentData
    lastParams::Vector{Float64}
    result::Vector{Float64}
    counter::Int
end


mutable struct CompiledExpression
    code::String
    funct::Function
    exprs::Vector{Expr}
    compInvolved::Vector{Symbol}
    counter::Int

    ldomain::AbstractLinearDomain
    domain::AbstractDomain
    cevals::HashVector{ComponentEvaluation}
end


mutable struct Model <: AbstractDict{Symbol, AbstractComponent}
    comp::HashVector{AbstractComponent}
    altValues::HashVector{Float64}
    compiled::Vector{CompiledExpression}
    results::Vector{Vector{Float64}}
    domainid::Vector{Int}
    Model() = new(HashVector{AbstractComponent}(), HashVector{Float64}(),
                  Vector{CompiledExpression}(), Vector{Vector{Float64}}(), Vector{Int}())
end


show(stream::IO, mime::MIME"text/plain", model::Model) = show(stream, model)
function show(stream::IO, model::Model)
    color = [229, 255]

    s = @sprintf "Components:\n"
    printstyled(color=:default, stream, s, bold=true)
    length(model) != 0  || (return nothing)

    s = @sprintf "%5s|%20s|%21s|%10s\n"  "#" "Component" "Type" "Alt. value"
    printstyled(color=:default, stream, s)

    color = sort(color)
    count = 0
    for (cname, comp) in model
        count += 1
        color = circshift(color, 1)

        ctype = split(string(typeof(comp)), ".")
        if ctype[1] == "ModelFit"
            ctype = ctype[2:end]
        end
        ctype = join(ctype[2:end], ".")

        altValue = model.altValues[cname]
        ss = (isfinite(altValue)  ?  @sprintf("%10.3g", altValue)  :  "")
        s = @sprintf "%5d|%20s|%21s|%10s\n" count string(cname) ctype ss
        printstyled(color=color[1], stream, s)
    end
    println(stream)

    printstyled(color=:default, stream, "Parameters:\n", bold=true)
    count = 0
    for (cname, comp) in model
        color = circshift(color, 1)
        count = show(stream, comp, cname=string(cname), count=count, color=color[1], header=false)
    end
    (length(model.compiled) == 0)  &&  (return)

    println(stream)
    s = @sprintf "Domains:\n"
    printstyled(color=:default, stream, s, bold=true)
    for i in 1:length(model.compiled)
        cc = model.compiled[i]

        s = @sprintf "#%d: " i
        printstyled(color=:default, stream, s, bold=true)
        show(stream, cc.domain)

        println(stream)
        s = @sprintf "%20s|%7s|%13s|%10s|%10s|%10s\n" "Component" "Counter" "Result size" "Min" "Max" "Mean"
        printstyled(color=:default, stream, s)

        color = sort(color)
        count = 0
        nonFinite = Vector{String}()

        for (cname, ceval) in cc.cevals
            count += 1
            color = circshift(color, 1)

            i = findall(isfinite.(ceval.result))
            v = view(ceval.result, i)

            s = @sprintf("%20s|%7d|%13s|%10.3g|%10.3g|%10.3g\n",
                         cname, ceval.counter, string(size(ceval.result)),
                         minimum(v), maximum(v), mean(v))
            printstyled(color=color[1], stream, s)

            nan = length(findall(isnan.(ceval.result)))
            inf = length(findall(isinf.(ceval.result)))
            if nan > 0  || inf > 0
                push!(nonFinite, @sprintf("%20s | NaN: %-10d   Inf: %-10d\n",
                                          cname, nan, inf))
            end
        end
        println(stream)
    end

    s = @sprintf "Expressions:\n"
    printstyled(color=:default, stream, s, bold=true)
    count = 0
    for i in 1:length(model.compiled)
        for j in 1:length(model.compiled[i].exprs)
            count += 1
            s = @sprintf "#%d (on domain #%d): %s\n" count i string(model.compiled[i].exprs[j])
            printstyled(color=:default, stream, s, bold=true)
        end
    end
    println(stream)

    s = @sprintf "%20s|%7s|%13s|%10s|%10s|%10s\n" "Expression" "Counter" "Result size" "Min" "Max" "Mean"
    printstyled(color=:default, stream, s)
    for i in 1:length(model.results)
        res = model.results[i]
        j = findall(isfinite.(res))
        v = view(res, j)

        s = @sprintf("%20s|%7d|%13s|%10.3g|%10.3g|%10.3g\n",
                     "#" * string(i), model.compiled[i].counter, string(size(res)),
                     minimum(v), maximum(v), mean(v))
        printstyled(color=:default, bold=true, stream, s)

        nonFinite = Vector{String}()
        nan = length(findall(isnan.(res)))
        inf = length(findall(isinf.(res)))
        if nan > 0  || inf > 0
            push!(nonFinite, @sprintf("%20s| NaN: %-10d   Inf: %-10d\n",
                                      "Meval", nan, inf))
        end

        if length(nonFinite) > 0
            println(stream)
            for s in nonFinite
                printstyled(color=:red, stream, bold=true, s)
            end
        end
    end
end


# ====================================================================
# FitParameter and FitResult structures, and associated show method.
#
struct FitParameter
    val::Float64
    unc::Float64
end


struct FitResult
    fitter::AbstractMinimizer
    param::HashVector{FitParameter}
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:Optimal, :NonOptimal, :Warn, :Error
    elapsed::Float64
end


function show(stream::IO, f::FitResult)
    #printstyled(color=:default, stream, "Minimizer:\n  " * string(typeof(f.fitter)), bold=true)
    #println(stream)
    #println(stream)

    color = [229, 255]

    s = @sprintf "Best fit values:\n"
    printstyled(color=:default, stream, s, bold=true)

    s = @sprintf "%5s|%20s|%10s|%10s|%10s|%10s\n"  "#" "Component" "Param." "Value" "Uncert." "Rel. unc. (%)"
    printstyled(color=:default, stream, s)

    color = sort(color)
    count = 0
    currentComp = ""
    for (pname, par) in f.param
        count += 1

        tmp = split(string(pname), "__")
        cname = tmp[1]
        pname = tmp[2]
        if currentComp != cname
            color = circshift(color, 1)
            currentComp = cname
        end

        s = @sprintf "%5d|%20s|%10s|%10.4g|%10.4g|%10.2g\n" count cname pname par.val par.unc par.unc/par.val*100.
        printstyled(color=color[1], stream, s)
    end

    println(stream)
    printstyled(color=:default, stream, "Summary:\n", bold=true)
    println(stream, @sprintf("  #Data  : %10d              Cost: %10.5g", f.ndata, f.cost))
    println(stream, @sprintf("  #Param : %10d              DOF : %10d", f.ndata-f.dof, f.dof))
    println(stream, @sprintf("  Elapsed: %10.4g s            Red.: %10.4g", f.elapsed, f.cost / f.dof))
    printstyled(color=:default, stream, "  Status :  ", bold=true)
    if f.status == :Optimal
        printstyled(color=:green, stream, "Optimal", bold=true)
    elseif f.status == :NonOptimal
        printstyled(color=:yellow, stream, "non-Optimal, see fitter output", bold=true)
    elseif f.status == :Warn
        printstyled(color=:yellow, stream, "Warning, see fitter output", bold=true)
    elseif f.status == :Error
        printstyled(color=:red, stream, "Error, see fitter output", bold=true)
    else
        printstyled(color=:magenta, stream, "Unknown (" * string(f.status) * "), see fitter output", bold=true)
    end
    println(stream)
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

compdata(domain::ModelFit.AbstractDomain, comp::FuncWrap) = FuncWrap_cdata(comp.func)
