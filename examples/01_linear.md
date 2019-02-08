# Example 1:

Fit a linear model with two parameters.


## Description

In this example we will generate a set of values from a linear model where we know the exact value of the model parameters (*true* values).  Then we will add some random noise to simulate a measurement process, and fit the obtained data with the linear model.  Finally, we will compare the resulting best fit values with the *true* parameter values.

Throughout the example we'll also describe the way data and model are summarized in the REPL by means of the `show` methods.

- **Components used**: `ScalarParam`
- **Functions used**: `Measures`, `Model`, `add_dom!`, `add_expr!`, `fit!`

To plot the results we will use the [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl) package, but any other will work as well.
```julia
using DataFitting, Gnuplot
```

## Data preparation
```julia
# "True" parameter values:
p0 = 1.
p1 = 2.

# Prepare the domain for the model (`x`) and calculate the "true"
# physical quantities (`y`):
x = 0.:0.1:10;
y = p0 .+ p1 .* x;

# Add some noise to simulate a measurement process
noise = 1.
y += noise .* randn(length(x));

# Wrap empirical measures and uncertainties in a `Measure` object
data = Measures(y, noise)
```
The last command will trigger the `show` method for the `Measure` object, and the following data will be shown in the REPL:
![show_data](https://github.com/gcalderone/DataFitting.jl/blob/master/examples/01_data.png)

All the `show` methods implemented in the DataFitting package are characterized by the blue frame around the tables (the appearence can be customized, see Example TODO).

In this case the first line tells you that that we are dealing with a `Measures_1D` (i.e. one-dimensional) object, with 101 samples.  This object has two properties: `val` and `unc`, each holding a vector of the value and uncertainties for each sample.  The columns in the table reports the minimu, maximum, average, median and standard deviation for the values and uncertainties vectors.

## Model preparation
Here we will show how to prepare the DataFitting model to fit the data.

```julia
# Prepare a model with two `ScalarParam` components, providing the guess values
model = Model(:a => ScalarParam(1), :b => ScalarParam(2))

# Set the domain to evaluate the model
add_dom!(model, x)

# Set the mathematical expression to evaluate the model
addexpr!(model, :(a .+ b .* domain[1]))

model
```
The last command will trigger the `show` method for the `Model` object, and the following data will be shown in the REPL:
![show_model](https://github.com/gcalderone/DataFitting.jl/blob/master/examples/01_model.png)

This shows the current status of the model.  From top to the bottom, it shows:
- The components used in the model, with their names (`a` and `b` in this case), type and a customary description.  The `F` column reports whether the component fitting is disabled (in this case all components are enabled).

- The parameters for all the components, with the components they belong to, the parameter names (`par` in this case), their current value, limits and customary description;

- The *"instrument(s)" used to obtain the empirical data.  *Instrument* here is a synonim for domain: there is one instrument for each domain added to the model with `add_dom`.  The user may add as many domains as needed.

- For each instrument the `show` method for a `Domain` object is called


# Fit the model to empirical data
result = fit!(model, data)

# Plot data (accessible through `data.val`), uncertainties
# (`data.unc`) and best fit model (`model[]`).  The independent
# variable is available through `model[].domain`.
using Gnuplot
@gp    model[].domain data.val data.unc "w yerr t 'Data'" :- 
@gp :- model[].domain model[].expr1 "w line t 'Model'"

# Compare "true" parameter values with best fit ones:
using Printf
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p0 result.a.par.val 100. * (result.a.par.val-p0) / p0
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p1 result.b.par.val 100. * (result.b.par.val-p1) / p1
