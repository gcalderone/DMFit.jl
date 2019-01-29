# ____________________________________________________________________
# Fall back methods
#
isfunction(comp::AbstractComponent) = false
description(comp::AbstractComponent) = ""
description(comp::AbstractComponent, param::Symbol) = ""
cdata(comp::AbstractComponent, domain::AbstractDomain) =
    error("Component " * string(typeof(comp)) * " must implement its own version of `cdata`.")
outsize(cdata::AbstractComponentData, domain::AbstractDomain) = length(domain)

# ____________________________________________________________________
# Built-in component: ScalarParam
#
struct ScalarParam <: AbstractComponent
    par::Parameter
    ScalarParam(val::Number) = new(Parameter(val))
end

struct ScalarParam_cdata <: AbstractComponentData; end
cdata(comp::ScalarParam, domain::AbstractDomain) = ScalarParam_cdata()

function evaluate!(cdata::ScalarParam_cdata, output::AbstractArray{Float64}, domain::AbstractDomain, par)
    output .= par
    return output
end


# ____________________________________________________________________
# Built-in component: FuncWrap
#
struct FuncWrap <: AbstractComponent
    func::Function
    p::Vector{Parameter}

    function FuncWrap(func::Function, args...)
        params = Vector{Parameter}()
        for i in 1:length(args)
            push!(params, Parameter(args[i]))
        end
        return new(func, params)
    end
end
description(c::FuncWrap) = "Function wrapper"

cdata(comp::FuncWrap, domain::AbstractDomain) = FuncWrap_cdata(comp.func)


# ____________________________________________________________________
# Built-in component: Smooth
#
struct Smooth <: AbstractComponent
    n::Int
    Smooth(n::Number) = new(Int(n))
end
isfunction(comp::Smooth) = true
description(c::Smooth) = "Smoothing function"

struct Smooth_cdata <: AbstractComponentData
    n::Int
    i::AbstractArray
end
cdata(comp::Smooth, domain::AbstractDomain) = Smooth_cdata(comp.n, 1:comp.n:(length(domain)-comp.n))
outsize(cdata::Smooth_cdata, domain::AbstractDomain) = length(cdata.i)

function evaluate!(cdata::Smooth_cdata, output::AbstractArray{Float64}, domain::AbstractDomain, v)
     output .= [mean(v[i:i+cdata.n-1]) for i in cdata.i]
    return output
end
