module DataFittingBasics

import DataFitting: AbstractDomain, Domain_1D, Domain_2D,
    Parameter, AbstractComponent, AbstractComponentData,
    cdata, evaluate!

include("OffsetSlope.jl")
include("Polynomial.jl")
include("Gaussian.jl")
include("Lorentzian.jl")
end
