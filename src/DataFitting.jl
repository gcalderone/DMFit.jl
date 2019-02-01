# ====================================================================
#                         MODULE DataFitting
# ====================================================================
module DataFitting

using Printf
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

export Domain, CartesianDomain, Measures,
    FuncWrap, ScalarParam, Smooth,
    Model, getparamvalues, setparamvalues!,
    addcomp!, setfixed!,
    add_dom!, rm_dom!, dom_count,
    addexpr!, replaceexpr!, setflag!,
    evaluate!, fit!, fit,
    test_component

const compsep = "_"

include("Types.jl")
include("components.jl")
include("private.jl")
include("userinterface.jl")
include("show.jl")

end
