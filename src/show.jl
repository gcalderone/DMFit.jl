function printmain(io::IO, args...)         ;                    printstyled(io, args..., "\n"; bold=true, color=:green); end
function printerr( io::IO, args...)         ;                    printstyled(io, args..., "\n"; bold=true, color=:red); end
function printtype(io::IO, args...)         ; print(io, "    "); printstyled(io, args..., "\n"; bold=true); end
#function printhead(io::IO, args...)         ; print(io, "   ╭"); printstyled(io, args..., "\n"; color=:underline, bold=true); end
#function printcell(io::IO, args...; u=false); print(io, "   ⋮"); printstyled(io, args..., "\n"; color=(u  ?  :underline  :  :default)); end


mutable struct MyTablePrinter
    head::String
    rule::String
    tail::String
    width::Int
    pendingrule::Bool
end
const mytable = MyTablePrinter("", "", "", 0, false)


function printhead(io::IO, args...)
    global mytable

    tmp = sprint(print, args...)
    ss = split(tmp, "│")
    for i in 1:length(ss)
        ss[i] = join(fill("─", length(ss[i])))
    end
    mytable.head = join(ss, "┬")
    mytable.rule = join(ss, "┼")
    mytable.tail = join(ss, "┴")
    mytable.width = length(tmp)+1
    mytable.pendingrule = false
    
    tmp = join(fill("─", mytable.width))
    printstyled(io, color=:light_blue, "    ╭", mytable.head, "╮\n    │")
    printstyled(io, bold=true, color=:light_blue, args..., "")
    printstyled(io, color=:light_blue, "│")
    println(io)
end
function printcell(io::IO, args...; lastingroup=false, tail=false)
    global mytable
    if mytable.pendingrule
        tmp = join(fill("─", mytable.width))
        printstyled(io, color=:light_blue, "    ├")
        print(io, mytable.rule)
        printstyled(io, color=:light_blue, "┤")
        println(io)
        mytable.pendingrule = false
    else
        mytable.pendingrule = lastingroup
    end
    printstyled(io, color=:light_blue, "    │")
    print(io, args...)
    tmp = length(sprint(print, args...))+1
    if tmp <= mytable.width
        tmp = join(fill(" ", mytable.width-tmp))
        printstyled(io, color=:light_blue, tmp, "│")
    end
    println(io)
end
function printtail(io::IO)
    global mytable
    tmp = join(fill("─", mytable.width))
    printstyled(io, color=:light_blue, "    ╰", mytable.tail, "╯\n")
    mytable.pendingrule = false
end


# │  ⃒  ◤
# ╭────┬─────╮
# │    │     │
# ├────┼─────┤
# ╰────┴─────╯


function show(io::IO, dom::AbstractCartesianDomain)
    printtype(io, typeof(dom), "  length: ", length(dom), "")
    s = @sprintf("%5s │ %8s │ %10s │ %10s │ %10s │ %10s",
                 "Dim.", "Size", "Min val", "Max val", "Min step", "Max step")
    printhead(io, s)
    for i in 1:ndims(dom)
        a = dom[i]
        b = 0
        if length(a) > 1
            b = a .- circshift(a, 1)
            b = b[2:end]
        end
        s = @sprintf("%5d │ %8d │ %10.4g │ %10.4g │ %10.4g │ %10.4g",
                     i, length(a),
                     minimum(a), maximum(a),
                     minimum(b), maximum(b))
        printcell(io, s)
    end
    printtail(io)
end


function show(io::IO, dom::AbstractLinearDomain)
    printtype(io, typeof(dom), " length: ", length(dom), "")
    s = @sprintf("%5s │ %10s │ %10s",
                 "Dim.", "Min val", "Max val")
    printhead(io, s)
    for i in 1:ndims(dom)
        s = @sprintf("%5d │ %10.4g │ %10.4g",
                     i, getaxismin(dom, i), getaxismax(dom, i))
        printcell(io, s)
    end
    printtail(io)
end


# Special case for Domain_1D: treat it as a Cartesian domain, despite it is a Linear one.
function show(io::IO, dom::Domain_1D)
    printtype(io, typeof(dom), " length: ", length(dom), "")
    s = @sprintf("%5s │ %8s │ %10s │ %10s │ %10s │ %10s",
                 "Dim.", "Size", "Min val", "Max val", "Min step", "Max step")
    printhead(io, s)
    a = dom[1]
    b = 0
    if length(a) > 1
        b = a .- circshift(a, 1)
        b = b[2:end]
    end
    s = @sprintf("%5d │ %8d │ %10.4g │ %10.4g │ %10.4g │ %10.4g",
                 1, length(a),
                 minimum(a), maximum(a),
                 minimum(b), maximum(b))
    printcell(io, s)
    printtail(io)
end


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
    printtail(io)
    
    if length(nonFinite) > 0
        println(io)
        for s in nonFinite
            printerr(s)
        end
    end
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

        printcell(io, s, lastingroup=(localcount == lastcount))
    end
    return count
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
    printtail(io)
    println(io)

    printmain(io, "Parameters:")
    count = 0
    for (cname, comp) in components(model)
        count = show(io, comp, cname=string(cname), count=count, header=false)
    end
    printtail(io)
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
            printcell(io, lastingroup=(localcount == lastcount),
                      @sprintf("%3d │ %10s │ %7d │ %10.3g │ %10.3g │ %10.3g │ %10d │ %10d │ %s",
                               ii, "Expr #"*string(jj), ce.counter,
                               minimum(v), maximum(v), mean(v), nan, inf, ce.exprs[jj]))
            countexpr += 1
        end
    end
    printtail(io)
    println(io)
    printmain(io, "Total expressions: ", countexpr)
end


function show(io::IO, comp::BestFitComp; count=0, cname="")
    if count == 0
        printhead(io, @sprintf "%5s │ %20s │ %10s │ %10s │ %10s │ %10s"  "#" "Component" "Param." "Value" "Uncert." "Rel.unc.(%)")
    end
    localcount = 0;  lastcount = length(getfield(comp, :params))
    for (pname, params) in getfield(comp, :params)
        localcount += 1
        if typeof(params) == Vector{BestFitParam}
            for ii in 1:length(params)
                count += 1
                par = params[ii]
                spname = string(pname) * "[" * string(ii) * "]"
                printcell(io, lastingroup=((localcount == lastcount)  &&  (ii == length(params))),
                                 @sprintf("%5d │ %20s │ %10s │ %10.4g │ %10.4g │ %10.2g", count, cname,
                                          spname, par.val, par.unc, par.unc/par.val*100.))
            end
        else
            count += 1
            par = params
            spname = string(pname)
            s = @sprintf("%5d │ %20s │ %10s │ %10.4g │ %10.4g │ %10.2g", count, cname,
                         spname, par.val, par.unc, par.unc/par.val*100.)
            printcell(io, s, lastingroup=(localcount == lastcount))
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
    printtail(io)

    println(io)
    printmain(io, "Summary:")
    println(io, @sprintf("    #Data  : %10d              Cost: %10.5g", f.ndata, f.cost))
    println(io, @sprintf("    #Param : %10d              DOF : %10d", f.ndata-f.dof, f.dof))
    println(io, @sprintf("    Elapsed: %10.4g s            Red.: %10.4g", f.elapsed, f.cost / f.dof))
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
