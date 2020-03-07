mutable struct Polynomial <: AbstractComponent
    coeff::Vector{Parameter}

    function Polynomial(args...)
        coeff = [Parameter(arg) for arg in args]
        return new(coeff)
    end
end

struct Polynomial_cdata <: AbstractComponentData; end
cdata(comp::Polynomial, domain::AbstractDomain) = Polynomial_cdata()

function evaluate!(cdata::Polynomial_cdata, output::AbstractArray{Float64}, domain::AbstractDomain, coeffs...)
    output .= coeffs[1]
    x = domain[1]
    for deg in 1:length(coeffs)-1
        output .+= x.^deg .* coeffs[deg+1]
    end
    return output
end
