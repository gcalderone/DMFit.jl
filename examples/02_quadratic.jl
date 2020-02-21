#= 
Example 2: Fit a linear model first, then add a further component and
fit a quadratic model.

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
y .+= noise .* randn(length(x));

# Wrap empirical measures and uncertainties in a `Measure` object
data = Measures(y, noise)

# Prepare the model.  Note that we start with an empty model, and add
# the components one by one:
model = Model()
add_comp!(model, :p0 => ScalarParam(1))    # add 1st component
add_comp!(model, :p1 => ScalarParam(2))    # add 2nd component
add_dom!(model, x)                         # set domain
add_expr!(model, :(p0 .+ p1 .* domain[1])) # add model expression

# Fit the model to empirical data
result = fit!(model, data)

# Plot data and best fit model
using Gnuplot
@gp    domain(model) data.val data.unc "w yerr t 'Data'" :- 
@gp :- domain(model) model() "w line t 'Model'"

# Print the reduced chi-squared and the probability that random noise
# may lead to higher chi-squared values.  Note that the latter is
# given as logarithmic value
println("Red. χ^2 = ", result.cost / result.dof)
println("Test probability = ", 10^result.log10testprob)

# Clearly this is not a good fit.  Let's add a 3rd `ScalarParam`
# component to obtain a quadratic model
add_comp!(model, :p2 => ScalarParam(3))          # add 3rd component

# Replace model expression
replaceexpr!(model, :(p0 .+ p1 .* domain[1] .+ p2 .* domain[1].^2))

# Fit model to empirical data and plot
result = fit!(model, data)
@gp    domain(model) data.val data.unc "w yerr t 'Data'" :- 
@gp :- domain(model) model() "w line t 'Model'"

# Compare "true" parameter values with best fit ones:
using Printf
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p0 result.p0.par.val 100. * (result.p0.par.val-p0) / p0
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p1 result.p1.par.val 100. * (result.p1.par.val-p1) / p1
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p2 result.p2.par.val 100. * (result.p2.par.val-p2) / p2
