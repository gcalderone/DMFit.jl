#= 
Example 3:  
  Fit the following analytical model with 3 parameters:
  (p1  +  p2 * x  +  p3 * x^2)  *  cos(x)

A `FuncWrap` component is used to wrap a common Julia function
=#
using Random, DataFitting

# Define the analytic model function
f(x, p1, p2, p3) = @. (p1  +  p2 * x  +  p3 * x^2)  *  cos(x)

# "True" parameter values:
params = [1, 2, 3];

# "True" physical quantities
x = 0:0.1:15
y = f(x, params...);

# Add some noise to simulate a measurement process
noise = 0.1 .* y
y .+= noise .* randn(length(x));

# Wrap empirical measures and uncertainties in a `Measure` object
data = Measures(y, noise)

# Prepare the model with a `FuncWrap` component.  Provide the guess
# values with the `params` vector.
model = Model(:comp1 => FuncWrap(f, params...))

# Set model domain and model expression, which in this case is simply
# the `comp1` component
add_dom!(model, x)
addexpr!(model, 1, :comp1)

# Fit the model to empirical data
result = fit!(model, data)

# Print best fit parameters
for i in 1:length(result.comp1.p)
    println("p[$i] = ", result.comp1.p[i].val, " ± ", result.comp1.p[i].unc)
end

# Plot data and best fit model
using Gnuplot
@gp    model(:domain) data.val data.unc "w yerr" :- 
@gp :- model(:domain) model() "w line"
