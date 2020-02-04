# Example 1:

Fit a linear model with two parameters.


## Description

In this example we will generate a set of values from a linear model where we know the exact value of the model parameters (*true* values).  Then we will add some random noise to simulate a measurement process, and fit the obtained data with the linear model.  Finally, we will compare the best fit values with the *true* parameter values.

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

In this case the first line reports the type of the object: `Measures_1D` (i.e. one-dimensional), with 101 samples.  This object has two properties: `val` and `unc`, each containing a vector for the value and uncertainties for each sample.  The columns in the table reports the minimum, maximum, average, median and standard deviation for the values and uncertainties vectors.

## Model preparation
```julia
# Prepare a model with two `ScalarParam` components, providing the guess values
model = Model(:a => ScalarParam(1), :b => ScalarParam(2))

# Set the domain to evaluate the model
add_dom!(model, x)

# Set the mathematical expression to evaluate the model
addexpr!(model, :(a .+ b .* domain[1]))

show(model)
```
The last command will trigger the `show` method for the `Model` object, and the following data will be shown in the REPL:
![show_model](https://github.com/gcalderone/DataFitting.jl/blob/master/examples/01_model.png)

This shows the current status of the model.  From top to bottom, it shows:
- The components used in the model, with their names (`a` and `b` in this case), type and a optional description.  The `F` column reports whether the component fitting is disabled. i.e. if the component is *fixed* (in this case all components are enabled).

- All the model parameters, with the components they belong to, their names (`par` in this case), their current value, limits and optional descriptions;

- The *"instrument(s)"* used to obtain the empirical data.  *Instrument* here is a synonim for domain: there is one instrument for each domain added to the model with `add_dom!`.  The user may add as many domains as needed.  For each domain the corresponding `show` method is invoked, and the domain details reported.  In this case the domain is one-dimensional, and has 101 samples.  For each dimension of the domain the size, minimum and maximum values, and minimum and maximum steps are reported;

- For each instrument (i.e., model domain) the list of evaluated expressions are reported. The entries above the horizontal white line are the components involved in the expressions (in this case `a` and `b`).  The entries below the lines are the expressions added with `addexpr!`.  For each expression the table reports how many times the expression has been evaluated, as well as the minimum, maximum and mean value of the results.  The last two columns reports respectively a flag to hgihlight NaN/infinite values in the evaluation, and the actual mathematical expression.  Note that the `expr1` expression has an arrow in the first column: it indicates that the results of this expression is to be compared with empirical data.  You may toggle this flag with `setflag!`;

- Finally, it shows the total number of instruments defined, and the number of datasets required for fitting.


## Fit the model to empirical data and produce plots
```julia
result = fit!(model, data)

# Plot data (accessible through `data.val`), uncertainties
# (`data.unc`) and best fit model (`model()`).  The independent
# variable is available through `model[].domain`.

@gp    model(:domain) data.val data.unc "w yerr t 'Data'" :- 
@gp :- model(:domain) model() "w line t 'Model'"

# Compare "true" parameter values with best fit ones:
using Printf
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p0 result.a.par.val 100. * (result.a.par.val-p0) / p0
@printf "%f ∼ %f  (accuracy: %5.3f%%)" p1 result.b.par.val 100. * (result.b.par.val-p1) / p1
```
![show_results](https://github.com/gcalderone/DataFitting.jl/blob/master/examples/01_results.png)
