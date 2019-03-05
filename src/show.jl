mutable struct ShowSettings
    tableformat::PrettyTableFormat
    floatformat::String
    border::Crayon
    header::Crayon
    subheader::Crayon
    fixed::Crayon
    error::Crayon
    section::Crayon
    fixedpars::Bool
    ShowSettings() = new(unicode_rounded, "%10.4g",
                         crayon"light_blue", crayon"light_blue bold",
                         crayon"dark_gray bold", crayon"dark_gray",
                         crayon"light_red blink", crayon"yellow bold",
                         true)
end

const showsettings = ShowSettings()

printtable(args...; kw...) = pretty_table(args..., showsettings.tableformat; kw...,
                                          border_crayon=showsettings.border,
                                          header_crayon=showsettings.header,
                                          subheader_crayon=showsettings.subheader)

section(io, args...) = println(io, showsettings.section, args...)

function left(s::String, maxlen::Int)
    (length(s) <= maxlen)  &&  (return s)
    return s[1:maxlen]
end

function show(io::IO, dom::AbstractCartesianDomain)
    section(io, "Cartesian domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    table = Matrix{Union{Int,Float64}}(undef, ndims(dom), 6)
    for i in 1:ndims(dom)
        a = dom[i]
        b = 0
        if length(a) > 1
            b = a .- circshift(a, 1)
            b = b[2:end]
        end
        table[i, 1] = i
        table[i, 2] = length(a)
        table[i, 3:6] = [minimum(a), maximum(a), minimum(b), maximum(b)]
    end
    printtable(io, table, ["Dim", "Size", "Min val", "Max val", "Min step", "Max step"],
               formatter=ft_printf(showsettings.floatformat, 3:6))
end


function show(io::IO, dom::AbstractLinearDomain)
    section(io, "Linear domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    table = Matrix{Union{Int,Float64}}(undef, ndims(dom), 3)
    for i in 1:ndims(dom)
        table[i, 1] = i
        table[i, 2] = getaxismin(dom, i)
        table[i, 3] = getaxismax(dom, i)
    end
    printtable(io, table, ["Dim", "Min val", "Max val"],
               formatter=ft_printf(showsettings.floatformat, 2:3))
end


# Special case for Domain_1D: treat it as a Cartesian domain, despite it is a Linear one.
function show(io::IO, dom::Domain_1D)
    section(io, "Domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    table = Matrix{Union{Int,Float64}}(undef, ndims(dom), 6)
    for i in 1:ndims(dom)
        a = dom[i]
        b = 0
        if length(a) > 1
            b = a .- circshift(a, 1)
            b = b[2:end]
        end
        table[i, 1] = i
        table[i, 2] = length(a)
        table[i, 3:6] = [minimum(a), maximum(a), minimum(b), maximum(b)]
    end
    printtable(io, table, ["Dim", "Size", "Min val", "Max val", "Min step", "Max step"],
               formatter=ft_printf(showsettings.floatformat, 3:6))
end


function show(io::IO, data::AbstractData)
    section(io, typeof(data), "   length: ", (length(data.val)))
    table = Matrix{Union{String,Float64}}(undef, 0, 7)

    names = fieldnames(typeof(data))
    for name in names
        a = getfield(data, name)
        nan = length(findall(isnan.(a))) + length(findall(isinf.(a)))
        a = a[findall(isfinite.(a))]
        table = vcat(table, [(nan > 0  ?  "⚠"  :  "") string(name) minimum(a) maximum(a) mean(a) median(a) std(a)])
    end
    printtable(io, table, ["⚠", "", "Min", "Max", "Mean", "Median", "Std. dev."],
               formatter=ft_printf(showsettings.floatformat, 3:7))
end


function preparetable(wcomp::WComponent)
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    comp = wcomp.comp
    cname = string(wcomp.cname)
    for (pname, param) in getparams(wcomp)
        (!showsettings.fixedpars)  &&  (wcomp.fixed  ||  param.fixed)  &&  continue
        ss = string(param._private.pname) * (param._private.index >= 1   ?   "["*string(param._private.index)*"]"  :  "")
        range = (param.fixed  ?  "FIXED"  :  strip(@sprintf("%7.2g:%-7.2g", param.low, param.high)))
        (range == "-Inf:Inf")  &&  (range = "")
        log = (param.log  ?  "LOG"  :  "")
        table = vcat(table, [cname ss param.val range log param.expr description(comp, param._private.pname)])
        push!(fixed, param.fixed)
        push!(error, !(param.low <= param.val <= param.high))
        cname = ""
    end
    return (table, fixed, error)
end

function preparetable(dict::OrderedDict{Symbol, WComponent})
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    hrule = Vector{Int}()
    for (cname, wcomp) in dict
        (t, f, e) = preparetable(wcomp)
        table = vcat(table, t)
        append!(fixed, f)
        append!(error, e)
        push!(hrule, length(error))
    end
    return (table, fixed, error, hrule)
end

show(io::IO, w::UI{WComponent}) = show(io, wrappee(w))
show(io::IO, wcomp::WComponent) = show(OrderedDict{Symbol, WComponent}(wcomp.cname => wcomp))

function show(io::IO, dict::OrderedDict{Symbol, WComponent})
    (table, fixed, error, hrule) = preparetable(dict)
    printtable(io, table , ["Component" "Param." "Value" "Range" "Log" "Expr" "Notes"], alignment=:l,
               hlines=hrule, formatter=ft_printf(showsettings.floatformat, [3]),
               highlighters=(Highlighter((data,i,j) -> fixed[i], showsettings.fixed),
                             Highlighter((data,i,j) -> (error[i] && j==5), showsettings.error)))
end

show(io::IO, mime::MIME"text/plain", model::Model) = show(io, model)
show(io::IO, w::UI{Model}) = show(io, wrappee(w))
function show(io::IO, model::Model)
    _evaluate!(model)
    section(io, "Model components:")
    length(model.comp) != 0  || (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 4)
    for (cname, wcomp) in model.comp
        ctype = split(string(typeof(wcomp.comp)), ".")
        (ctype[1] == "DataFitting")  &&   (ctype = ctype[2:end])
        ctype = join(ctype, ".")
        table = vcat(table, [string(cname) (wcomp.fixed  ?  "F"  :  "") ctype description(wcomp.comp)])
    end
    printtable(io, table, ["Component" "F" "Type" "Description"], alignment=:l)

    println(io)
    section(io, "Parameters:")
    show(io, model.comp)

    if length(model.instruments) == 0
        println(io)
        section(io, "Instrument(s): 0")
        return nothing
    end

    for i in 1:length(model.instruments)
        println(io)
        section(io, "Instrument #$i ")
        show(io, model.instruments[i])
    end

     tmp = length(model.index1d) - 1
     (tmp < 0)  &&  (tmp = 0)
     section(io, "Instrument(s): ", length(model.instruments),
               ".  Dataset(s) required for fitting: ", tmp)
end



function show(io::IO, instr::Instrument)
#     if ((isa(io, Base.TTY)  ||  isa(io, IOContext))  &&  (displaysize(io)[2] >= 80))
    show(io, instr.domain)
    println(io)

    (length(instr.compevals) == 0)  &&  (return nothing)
    section(io, "Evaluated expression(s):")
    table = Matrix{Union{String,Int,Float64}}(undef, length(instr.compevals)+length(instr.exprs), 8)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    hrule = Vector{Int}()

    for i in 1:length(instr.compevals)
        cname = instr.compnames[i]
        ceval = instr.compevals[i]
        result = ceval.result
        v = view(result, findall(isfinite.(result)))
        (length(v) == 0)  &&  (v = [NaN])
        nan = length(findall(isnan.(result)))
        inf = length(findall(isinf.(result)))
        table[i, 1] = ""
        table[i, 2] = string(cname)
        table[i, 3] = ceval.counter
        table[i, 4:6] = [minimum(v), maximum(v), mean(v)]
        table[i, 7] = (nan+inf > 0 ? "⚠" : "")
        table[i, 8] = ""
        push!(fixed, ceval.fixed)
        push!(error, (nan+inf > 0))
    end
    i0 = length(instr.compevals)
    
    for i in 1:length(instr.exprs)
        result = instr.exprevals[i]
        v = view(result, findall(isfinite.(result)))
        nan = length(findall(isnan.(result)))
        inf = length(findall(isinf.(result)))
        table[i0+i, 1] = (instr.exprcmp[i]  ?  "⇒"  :  "")
        table[i0+i, 2] = string(instr.exprnames[i])
        table[i0+i, 3] = instr.counter
        table[i0+i, 4:6] = [minimum(v), maximum(v), mean(v)]
        table[i0+i, 7] = (nan+inf > 0 ? "⚠" : "")
        table[i0+i, 8] = left(string(instr.exprs[i]), 30)
        push!(fixed, false)
        push!(error, (nan+inf > 0))
    end

    printtable(io, table, ["", "Label", "Counter", "Min", "Max", "Mean", "⚠", "Expr"], alignment=:l,
               hlines=[i0], formatter=ft_printf(showsettings.floatformat, 4:6),
               highlighters=(Highlighter((data,i,j) -> fixed[i], showsettings.fixed),
                             Highlighter((data,i,j) -> (error[i] && j==5), showsettings.error)))
    println(io)
end


function show(io::IO, par::Parameter)
    if par.fixed
        println(io, "Value : ", par.val, "   (FIXED)")
    else
        println(io, "Value : ", par.val, "  [", par.low , " : ", par.high, "]")
        if par.expr != ""
            println(io, "Expr : ", par.expr)
        end
    end
    if par.log
        println(io, "(use logarithmic value in fit)")
    end
end

show(io::IO, par::FitParam) = println(io, par.val, " ± ", par.unc)

show(io::IO, w::UI{FitComp}) = show(io, wrappee(w))


function preparetable(comp::FitComp)
    table = Matrix{Union{String,Float64}}(undef, 0, 4)
    fixed = Vector{Bool}()
    error = Vector{Bool}()

    for (pname, params) in comp.params
        if typeof(params) == Vector{FitParam}
            for ii in 1:length(params)
                par = params[ii]
                (!showsettings.fixedpars)  &&  (par.fixed)  &&  continue
                spname = string(pname) * "[" * string(ii) * "]"
                table = vcat(table, ["" spname par.val par.unc])
                push!(fixed, par.fixed)
                push!(error, !isfinite(par.unc))
            end
        else
            par = params
            (!showsettings.fixedpars)  &&  (par.fixed)  &&  continue
            spname = string(pname)
            table = vcat(table, ["" spname par.val par.unc])
            push!(fixed, par.fixed)
            push!(error, !isfinite(par.unc))
        end
    end
    return (table, fixed, error)
end


function show(io::IO, w::UI{FitResult})
    res = wrappee(w)
    section(io, "Best Fit results:")

    table = Matrix{Union{String,Float64}}(undef, 0, 4)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    hrule = Vector{Int}()    
    for (cname, comp) in res.bestfit
        if length(comp.params) > 0
            (t, f, e) = preparetable(comp)
            t[1,1] = string(cname)
            table = vcat(table, t)
            append!(fixed, f)
            append!(error, e)
            push!(hrule, length(error))
        end
    end
    printtable(io, table , ["Component" "Param." "Value" "Uncert."], alignment=:l,
               hlines=hrule, formatter=ft_printf(showsettings.floatformat, [3,4]),
               highlighters=(Highlighter((data,i,j) -> fixed[i], showsettings.fixed),
                             Highlighter((data,i,j) -> (error[i]  &&  (!fixed[i])  &&  (j==4)), showsettings.error)))
    
    println(io)
    println(io, @sprintf("    #Data  : %10d              Cost : %-10.5g", res.ndata, res.cost))
    println(io, @sprintf("    #Param : %10d              Red. : %-10.4g", res.ndata-res.dof, res.cost / res.dof))
    print(  io, @sprintf("    DOF    : %10d              ", res.dof))
    if res.log10testprob < -3
        println(io, @sprintf("Prob.: 10^%-10.4g", res.log10testprob))
    else
        println(io, @sprintf("Prob.: %10.4g", 10^res.log10testprob))
    end
    printstyled(io, "    Status :  ")
    if res.status == :Optimal
        printstyled(color=:green, io, @sprintf("%-15s", "Optimal"))
    elseif res.status == :NonOptimal
        printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Non Optimal"))
    elseif res.status == :Warn
         printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Warning"))
    elseif res.status == :Error
        printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Error"))
    else
        printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Unknown (" * string(res.status) * "), see fitter output"))
    end
    println(io, @sprintf("        Elapsed: %-10.4g s", res.elapsed))
end

