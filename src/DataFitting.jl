# ====================================================================
#                         MODULE DataFitting
# ====================================================================
module DataFitting

using Printf
using Statistics
using DataStructures

import Base.show
import Base.ndims
import Base.size
import Base.length
import Base.getindex
import Base.propertynames
import Base.getproperty
import Base.append!

export Domain, CartesianDomain, Measures,
    FuncWrap, ScalarParam, Smooth,
    Model,
    addcomp!, setenabled!,
    add_dom!, rm_dom!, dom_count, recompile!,
    addexpr!, replaceexpr!, setflag!,
    evaluate!, fit!,
    test_component

const compsep = "_"

include("Types.jl")
include("components.jl")
@code_ndim 1
@code_ndim 2

include("private.jl")
include("userinterface.jl")
include("show.jl")

end
