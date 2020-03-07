module DataFitting

using Printf, PrettyTables
using Statistics, Distributions
using DataStructures

import Base.show
import Base.ndims
import Base.size
import Base.length
import Base.getindex
import Base.propertynames
import Base.getproperty
import Base.setproperty!
import Base.append!
import Base.isequal
import Base.reshape

export Domain, CartesianDomain, Measures,
    FuncWrap, ScalarParam, Smooth,
    Model, getparamvalues, setparamvalues!, resetcounters!, componenttype, domain,
    add_comp!, add_dom!, rm_dom!, dom_count, add_expr!, replaceexpr!, setflag!,
    evaluate!, fit!, probe,
    test_component

const compsep = "_"

include("Types.jl")
include("components.jl")
include("private.jl")
include("userinterface.jl")
include("show.jl")

const extfunc = Dict{Symbol, Function}()

end
