#= 
Example 2:  
  Fit the following analytical model with 5 parameters:
  (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)

A `FuncWrap` component is used to wrap a common Julia function
=#
using Random, DataFitting

# Define the analytic model function
f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)

# "True" parameter values:
params = [1, 1.e-3, 1.e-6, 4, 5]

# "True" values 
x = 1.:5:5000
y = f(x, params...);

# Add some noise to simulate a measurement process
noise = 1.
y .+= noise .* randn(MersenneTwister(0), length(x));

# Prepare:
# - domain for model evaluation
# - empirical data
data = Measures(y, noise)

# Build a model with appropriate components and parameter guess values
model = Model(:comp1 => FuncWrap(f, params...))

# Add an instrument (on the model domain), and an expression involving
# defined components:
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
@gp    model[].domain data.val :- 
@gp :- model[].domain model[].expr1 "w line"
