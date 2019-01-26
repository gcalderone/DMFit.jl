#= 
Example 1:
Fit a linear model.

A linear model is prepared using two `ScalarParam` components, and fit
to empirical data.
=#
using Random, DataFitting

# "True" parameter values:
p0 = 1.
p1 = 2.

# Evaluate the model with "True" parameters
x = 0.:0.1:10;
y = collect(p0 .+ p1 .* x);

# Add some noise to simulate a measurement process
noise = 1.
y .+= noise .* randn(MersenneTwister(0), length(x));

# Wrap empirical measure and uncertainties in a `Measure` object
data = Measures(y, noise)

# Prepare a model with two `ScalarParam` components, add a model
# domain and an expression involving components:
model = Model(:p0 => ScalarParam(p0), :p1 => ScalarParam(p1))
add_dom!(model, x)                        # add domain
addexpr!(model, :(p0 .+ p1 .* domain[1])) # add model expression

# Fit the model to empirical data
result = fit!(model, data)

# Plot data and best fit model
using Gnuplot
@gp    model[].domain data.val data.unc "w yerr t 'Data'" :- 
@gp :- model[].domain model[].expr1 "w line t 'Model'"

# Compare "true" parameter values with best fit ones:
using Printf
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p0 result.p0.par.val 100. * (result.p0.par.val-p0) / p0
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p1 result.p1.par.val 100. * (result.p1.par.val-p1) / p1



# ====================================================================
# Components can be added directly in Model constructor, e.g.

