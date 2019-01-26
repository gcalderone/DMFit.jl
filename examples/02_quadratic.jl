#= 
Example 2: 
Fit a linear model first, then add a further component and fit a quadratic model.

This example is very similar to #1, but it shows how to check for
goodness of fit, and how to add a new component to obtain a better fit.
=#
using Random, DataFitting

# "True" parameter values:
p0 = 1.
p1 = 2.
p2 = 3.

# Evaluate the model with "True" parameters
x = 0.:0.1:10;
y = collect(p0 .+ p1 .* x .+ p2 .* x.^2);

# Add some noise to simulate a measurement process
noise = 1.
y .+= noise .* randn(MersenneTwister(0), length(x));

# Wrap empirical measure and uncertainties in a `Measure` object
data = Measures(y, noise)

# Prepare a model with just two `ScalarParam` components, add a model
# domain and an expression involving components:
model = Model()                           # Start with empty Model
addcomp!(model, :p0 => ScalarParam(p0))   # add 1st component
addcomp!(model, :p1 => ScalarParam(p1))   # add 2nd component
add_dom!(model, x)                        # add domain
addexpr!(model, :(p0 .+ p1 .* domain[1])) # add model expression

# Fit the model to empirical data
result = fit!(model, data)

# Plot data and best fit model
using Gnuplot
@gp    model[].domain data.val data.unc "w yerr t 'Data'" :- 
@gp :- model[].domain model[].expr1 "w line t 'Model'"

# Print the reduced chi-squared
println("Red. χ^2 = ", result.cost / result.dof)

# Clearly this is not a good fit.  Let's add a 3rd `ScalarParam`
# component to obtain a quadratic model
addcomp!(model, :p2 => ScalarParam(p2))   # add 3rd component
addexpr!(model, :(expr1 .+ p2 .* domain[1].^2)) # add model expression
setflag!(model, :expr1, false)

result = fit!(model, data)
@gp    model[].domain data.val data.unc "w yerr t 'Data'" :- 
@gp :- model[].domain model[].expr2 "w line t 'Model'"


# Compare "true" parameter values with best fit ones:
using Printf
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p0 result.p0.par.val 100. * (result.p0.par.val-p0) / p0
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p1 result.p1.par.val 100. * (result.p1.par.val-p1) / p1
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p2 result.p2.par.val 100. * (result.p2.par.val-p2) / p2






model = Model(:p0 => ScalarParam(p0), :p1 => ScalarParam(p1))
add_dom!(model, x)                        # add domain
addexpr!(model, :(p0 .+ p1 .* domain[1])) # add model expression
addcomp!(model, :p2 => ScalarParam(p2))   # add 3rd component
replaceexpr!(model, :expr1, :(p0 .+ p1 .* domain[1] .+ p2 .* domain[1].^2))
result = fit!(model, data)
@gp    model[].domain data.val data.unc "w yerr t 'Data'" :- 
@gp :- model[].domain model[].expr1 "w line t 'Model'"
