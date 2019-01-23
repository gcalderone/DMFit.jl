mutable struct PrintSettings
    prefix::String
    head::String
    rule::String
    tail::String
    width::Int
    pendingrule::Bool
    compact::Bool
    colormain::Symbol
    colortable::Symbol
end
const ps = PrintSettings("", "", "", "", 0, false, false, :magenta, :light_blue)

function showcompact(b::Bool=true)
    global ps
    ps.compact = b
end

function printbold()
    global ps
    (ps.compact)  &&  return false
    return true
end
function printcolormain()
    global ps
    (ps.compact)  &&  return :default
    return ps.colormain
end
function printcolortable()
    global ps
    (ps.compact)  &&  return :default
    return ps.colortable
end
function printerr( io::IO, args...)
    printstyled(io, args..., "\n"; bold=printbold(), color=:red)
end
function printmain(io::IO, args...; newline=true)
    printstyled(io, args...; bold=printbold(), color=printcolormain())
    (newline)  &&  (println(io))
end
function printsub(io::IO, args...; newline=true)
    global ps
    print(io, ps.prefix)
    printstyled(io, args...; bold=printbold())
    (newline)  &&  (println(io))
end

function printhead(io::IO, args...)
    global ps
    (ps.compact)  &&  (return println(io, args...))
    tmp = sprint(print, args...)
    ss = split(tmp, "│")
    for i in 1:length(ss)
        ss[i] = join(fill("─", length(ss[i])))
    end
    ps.head = join(ss, "┬")
    ps.rule = join(ss, "┼")
    ps.tail = join(ss, "┴")
    ps.width = length(tmp)+1
    ps.pendingrule = false

    tmp = join(fill("─", ps.width))
    printstyled(io, color=printcolortable(), ps.prefix, "╭", ps.head, "╮\n", ps.prefix, "│")
    printstyled(io, bold=printbold(), color=printcolortable(), args...)
    printstyled(io, color=printcolortable(), "│")
    println(io)
end
function printcell(io::IO, args...; lastingroup=false, tail=false)
    global ps
    (ps.compact)  &&  (return println(io, args...))
    if ps.pendingrule
        tmp = join(fill("─", ps.width))
        printstyled(io, color=printcolortable(), ps.prefix, "├")
        print(io, ps.rule)
        printstyled(io, color=printcolortable(), "┤")
        println(io)
        ps.pendingrule = false
    else
        ps.pendingrule = lastingroup
    end
    printstyled(io, color=printcolortable(), ps.prefix, "│")
    printstyled(io, args...)#, color=(lastingroup ? :underline : :default))
    tmp = length(sprint(print, args...))+1
    if tmp <= ps.width
        tmp = join(fill(" ", ps.width-tmp))
        printstyled(io, color=printcolortable(), tmp, "│")
    end
    println(io)
end
function printtail(io::IO)
    global ps
    (ps.compact)  &&  (return nothing)
    tmp = join(fill("─", ps.width))
    printstyled(io, color=printcolortable(), ps.prefix, "╰", ps.tail, "╯\n")
    ps.pendingrule = false
end

function show(io::IO, dom::AbstractCartesianDomain)
    printsub(io, "Cartesian domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    s = @sprintf("%3s │ %8s │ %10s │ %10s │ %10s │ %10s",
                 "Dim", "Size", "Min val", "Max val", "Min step", "Max step")
    printhead(io, s)
    for i in 1:ndims(dom)
        a = dom[i]
        b = 0
        if length(a) > 1
            b = a .- circshift(a, 1)
            b = b[2:end]
        end
        s = @sprintf("%3d │ %8d │ %10.4g │ %10.4g │ %10.4g │ %10.4g",
                     i, length(a),
                     minimum(a), maximum(a),
                     minimum(b), maximum(b))
        printcell(io, s)
    end
    printtail(io)
end


function show(io::IO, dom::AbstractLinearDomain)
    printsub(io, "Linear domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    s = @sprintf("%3s │ %10s │ %10s",
                 "Dim", "Min val", "Max val")
    printhead(io, s)
    for i in 1:ndims(dom)
        s = @sprintf("%3d │ %10.4g │ %10.4g",
                     i, getaxismin(dom, i), getaxismax(dom, i))
        printcell(io, s)
    end
    printtail(io)
end


# Special case for Domain_1D: treat it as a Cartesian domain, despite it is a Linear one.
function show(io::IO, dom::Domain_1D)
    printsub(io, "Domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    s = @sprintf("%3s │ %8s │ %10s │ %10s │ %10s │ %10s",
                 "Dim", "Size", "Min val", "Max val", "Min step", "Max step")
    printhead(io, s)
    a = dom[1]
    b = 0
    if length(a) > 1
        b = a .- circshift(a, 1)
        b = b[2:end]
    end
    s = @sprintf("%3d │ %8d │ %10.4g │ %10.4g │ %10.4g │ %10.4g",
                 1, length(a),
                 minimum(a), maximum(a),
                 minimum(b), maximum(b))
    printcell(io, s)
    printtail(io)
end


function show(io::IO, data::AbstractData)
    printsub(io, typeof(data), "   length: ", (length(data.measure)))
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
    printtail(io)

    if length(nonFinite) > 0
        println(io)
        for s in nonFinite
            printerr(s)
        end
    end
end


function show(io::IO, comp::AbstractComponent; header=true, count=0, cname="")
    (header)  &&  (printsub(io, typeof(comp)))
    if count == 0
        printhead(io, @sprintf "%-15s │ %-10s │ %10s │ %10s │ %10s │ %s"  "Component" "Param." "Value" "Low" "High" "Notes")
    end

    localcount = 0; lastcount = length(getparams(comp))
    for (pname, wparam) in getparams(comp)
        localcount += 1
        par = wparam.par
        note = ""
        (par.fixed)  &&  (note *= "FIXED")
        (par.expr != "")  &&  (note *= " expr=" * par.expr)
        count += 1
        s = @sprintf("%-15s │ %-10s │ %10.3g │ %10.3g │ %10.3g │ %s",
                     cname,
                     (wparam.index >= 1  ?  Symbol(wparam.pname, "[", wparam.index, "]")  :  pname),
                     par.val, par.low, par.high, note)

        printcell(io, s, lastingroup=(localcount == lastcount))
    end
    return count
end


show(io::IO, mime::MIME"text/plain", model::Model) = show(io, model)
show(io::IO, w::Wrap{Model}) = show(io, wrappee(w))
function show(io::IO, model::Model)
    printmain(io, "Model components:")
    length(model.comp) != 0  || (return nothing)

    printhead(io, @sprintf "%-15s │ %-30s"  "Component" "Type")
    count = 0
    for (cname, comp) in model.comp
        count += 1
        ctype = split(string(typeof(comp)), ".")
        (ctype[1] == "DataFitting")  &&   (ctype = ctype[2:end])
        ctype = join(ctype, ".")

        s = @sprintf "%-15s │ %-30s" string(cname) ctype
        printcell(io, s)
    end
    printtail(io)
    println(io)

    printsub(io, "Parameters:")
    count = 0
    for (cname, comp) in model.comp
        count = show(io, comp, cname=string(cname), count=count, header=false)
    end
    printtail(io)
    println(io)

    if length(model.instruments) == 0
        printmain(io, "Total instruments: 0")
        return nothing
    end

    printmain(io, "Instrument(s) defined: ", length(model.instruments))
    countcmp = 0
    for ii in 1:length(model.instruments)
        instr  = model.instruments[ii]
        countcmp += length(findall(instr.exprcmp))
    end
    printmain(io, "Dataset(s) required to fit this model: ", countcmp)
end

show(io::IO, w::Wrap{Instrument}) = show(io, wrappee(w))


function show(io::IO, instr::Instrument)
    function left(s::AbstractString, maxlen::Int)
        if maxlen < length(s)
            return s[1:maxlen-1] * "…"
        end
        return s
    end

    # Check max length of expressions
    exprmaxlength = 0
    for e in instr.exprs
        l = length(string(e))
        (exprmaxlength < l)  &&  (exprmaxlength = l)
    end
    #Check available space to fill the terminal
    availlength = exprmaxlength
    if (isa(io, IOContext))  &&  (displaysize(io)[2] >= 80)
        availlength = displaysize(io)[2]-76-length(ps.prefix)
    end
    (availlength > exprmaxlength)  &&  (availlength = exprmaxlength)
    (availlength < 4)  &&  (availlength = 4) # Leave space for "Expr" header

    printmain(io, "Instrument (", instr.label, ") on ", newline=false)
    show(io, instr.domain)
    println(io)

    printsub(io, "Evaluated expression(s):")
    printhead(io, @sprintf "%2s%-15s │ %7s │ %10s │ %10s │ %10s │ %1s │ %s " "" "Label" "Counter" "Min" "Max" "Mean" "⚠" "Expr"*join(fill(" ", availlength-4)))

    for jj in 1:length(instr.compevals)
        cname = instr.compnames[jj]
        ceval = instr.compevals[jj]

        result = ceval.result
        v = view(result, findall(isfinite.(result)))
        nan = length(findall(isnan.(result)))
        inf = length(findall(isinf.(result)))
        printcell(io, @sprintf("%2s%-15s │ %7d │ %10.3g │ %10.3g │ %10.3g │ %1s │ ",
                               "", cname, ceval.counter,
                               minimum(v), maximum(v), mean(v), (nan+inf > 0 ? "⚠" : "")),
                  lastingroup=(jj==length(instr.compevals)))
    end

    localcount = 0; lastcount = length(instr.exprs)
    for jj in 1:length(instr.exprs)
        localcount += 1
        result = instr.exprevals[jj]
        v = view(result, findall(isfinite.(result)))
        nan = length(findall(isnan.(result)))
        inf = length(findall(isinf.(result)))
        printcell(io, lastingroup=(localcount == lastcount),
                  @sprintf("%-2s%-15s │ %7d │ %10.3g │ %10.3g │ %10.3g │ %1s │ %s",
                           (instr.exprcmp[jj]  ?  "⇒"  :  ""),
                           instr.exprnames[jj], instr.counter,
                           minimum(v), maximum(v), mean(v), (nan+inf > 0 ? "⚠" : ""),
                           left(string(instr.exprs[jj]), availlength)))
    end
    printtail(io)
    println(io)
end


show(io::IO, w::Wrap{FitComp}) = show(io, wrappee(w))
function show(io::IO, comp::FitComp; header=true, cname="")
    if header
        printhead(io, @sprintf "%-15s │ %-10s │ %10s │ %10s │ %10s"  "Component" "Param." "Value" "Uncert." "Rel.unc.(%)")
    end
    localcount = 0;  lastcount = length(comp.params)
    for (pname, params) in comp.params
        localcount += 1
        if typeof(params) == Vector{FitParam}
            for ii in 1:length(params)
                par = params[ii]
                spname = string(pname) * "[" * string(ii) * "]"
                printcell(io, lastingroup=((localcount == lastcount)  &&  (ii == length(params))),
                                 @sprintf("%-15s │ %-10s │ %10.4g │ %10.4g │ %10.2g", cname,
                                          spname, par.val, par.unc, par.unc/par.val*100.))
            end
        else
            par = params
            spname = string(pname)
            s = @sprintf("%-15s │ %-10s │ %10.4g │ %10.4g │ %10.2g", cname,
                         spname, par.val, par.unc, par.unc/par.val*100.)
            printcell(io, s, lastingroup=(localcount == lastcount))
        end
    end
end


function show(io::IO, w::Wrap{FitResult})
    res = wrappee(w)
    printmain(io, "Best Fit results:")

    first = true
    for (cname, comp) in res.bestfit
        show(io, comp, header=first, cname=string(cname))
        first = false
    end
    printtail(io)

    println(io)
    println(io, @sprintf("    #Data  : %10d              Cost: %10.5g", res.ndata, res.cost))
    println(io, @sprintf("    #Param : %10d              DOF : %10d", res.ndata-res.dof, res.dof))
    println(io, @sprintf("    Elapsed: %10.4g s            Red.: %10.4g", res.elapsed, res.cost / res.dof))
    printstyled(io, "    Status :  ", bold=printbold())
    if res.status == :Optimal
        printstyled(color=:green, io, "Optimal", bold=printbold())
    elseif res.status == :NonOptimal
        printstyled(color=:yellow, io, "non-Optimal, see fitter output", bold=printbold())
    elseif res.status == :Warn
        printstyled(color=:yellow, io, "Warning, see fitter output", bold=printbold())
    elseif res.status == :Error
        printstyled(color=:red, io, "Error, see fitter output", bold=printbold())
    else
        printstyled(color=:magenta, io, "Unknown (" * string(res.status) * "), see fitter output", bold=printbold())
    end
    println(io)
end
