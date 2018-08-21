# ModelFit

# DMFit
Data-driven Model Fitting.

[![Build Status](https://travis-ci.org/gcalderone/ModelFit.jl.svg?branch=master)](https://travis-ci.org/gcalderone/ModelFit.jl)

`ModelFit` is a general purpose modeling and fitting framework for Julia.

**The package documentation is being updated...**


## Key features:
- it handles data of any dimensionality;
- the fitting model is evaluated on a user provided (Cartesian or Linear) domain;
- the fitting model is built up by individual *components*, either provided in the companion package by `ModelFit_Components` or implemented by the user.  All components are combined to evaluate the final model with a standard Julia mathematical expression;
- all components results are cached, so that repeated evaluations with the same parameters do not involve further calculations.  This is very important to speed up the fitting process when many components are involved;
- it allows fitting multiple data sets;
- it easily allows to use different minimizers, and compare their results and performances.  Currently two minimizers are supported ([LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) and [CMPFit](https://github.com/gcalderone/CMPFit.jl));
- provides several facilities for interactive fitting and results displaying.

## Installation:
On Julia v0.6:
```julia
Pkg.clone("https://github.com/gcalderone/ModelFit.jl.git")
```

On Julia v0.7/v1.0:
```julia
Pkg.develop("https://github.com/gcalderone/ModelFit.jl.git")
```

## Simple example:
Assume the model to be compared with empirical data has 5 parameters and the foolowing analytical formula:
```julia
f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)
```

To simulate a measurement process we'll evaluate the model on a domain and add some random noise to a model realization:

```julia
# Actual model parameters:
params = [1, 1.e-3, 1.e-6, 4, 5]

# Domain for model evaluation
x = 1:0.05:10000

# Evaluated model
y = f(x, params...);

# Random noise
rng = MersenneTwister(0);
noise = randn(rng, length(x));
```

In order to use the `ModelFit` framework we must create a `Domain` and a `Measures` objects to encapsulate the model domain and empirical data:
```julia
using ModelFit
dom = Domain(x)
data = Measures(y + noise, 1.)
```
The second argument to the `Measures` function is the uncertainty associated to each data sample.  It can either be an array with the same shape as the first argument or a scalar.  In the latter case all data samples will have the same uncertainty.

Also, we must create a `Model` object containing a reference to the analytical formula, and prepare it for evaluation on the domain `dom`
```julia
model1 = Model(:comp1 => FuncWrap(f, params...))
prepare!(model1, dom, :comp1)
```
The `comp1` symbol is the name we chose to identify the component in the model.  Any other valid symbol could have been used.

The parameter initial values are those given in the component constructors.  Such values can be changed as follows:
```julia
 model1[:comp1].p1.val = 1
 model1[:comp1].p2.val = 1.e-3
 ```
 etc.


Finally, we are ready to fit the model against empirical data:
```julia
result1 = fit!(model1, data)
```
Note that all the objects created so far (i.e. `dom`, `data`, `model1` and `result1`) have a dedicated `show` method to quickly and easily inspect their contents.

The procedure outlined above may seem cumbersome at first, however it is key to define very complex models and to improve performances, as shown below.



## Multiple components

A model is typically made up of many *components*, joined toghether with a standard mathematical expression.  The previous example can easily be re-implemented by splitting the analytical formula in 3 parts:
```julia
f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @. p4 * sin(p5 * x)
f3(x) = cos.(x)

model2 = Model(:comp1 => FuncWrap(f1, p1, p2, p3),
               :comp2 => FuncWrap(f2, p4, p5),
               :comp3 => FuncWrap(f3))
prepare!(model2, dom, :((comp1 + comp2) * comp3))
result2 = fit!(model2, data)
```
Now we used 3 components (named `comp1`, `comp2` and `comp3`) and combined them with the mathematical expression in the last argument of the call to `prepare!`. **Any** valid mathematical expression is allowed.

Note that we obtained exactly the same result as before, but we **gained a factor ~2** in execution time.  Such perfromance improvement is due to the fact that the component evaluations are internally cached, and are re-evaluated only if necessary, i.e. when one of the parameter value is modified by the minimzer algorithm. In this particular case the `comp3` component, having no free parameter, is evaluated only once.

To check how many time each component and the model are evaluated simply type the name of the `Model` object in the REPL or call the `show` method, i.e.: `show(model2)`, and look at the `Counter` column.


## Fitting multiple data sets
`ModelFit` allows to fit multiple data sets simultaneously.

Suppose you observe the same phenomenon with two different instruments and wish to use both data sets to constrain the model parameters.  Here we will simulate a second data sets by adding random noise to the previously calculated model, and creating a second `Measures` object: 
```julia
noise = randn(rng, length(x));
data2 = Measures(1.3 * (y + noise), 1.3)
```
Note that we multiplied all data by a factor 1.3, to simulate a different calibration between the instruments.  To take into account such calibration we push a further component into the model, named `calib`:
```julia
push!(model2, :calib, SimpleParam(1))
```
and prepare the model to evaluate a second expression:
```julia
prepare!(model2, dom, :(calib * ((comp1 + comp2) * comp3)))
```
Note that a call to `prepare!` modifies the `Model` object by adding a new (possibly different) expression to be evaluated.  If needed, also the `Domain` object used for the second expression may be different from the first one.  The latter posssibility allows, for instance, to fit multiple data sets each observed with different instruments.

Now we can fit both data sets as follows:
```julia
result2 = fit!(model2, [data, data2])
```

## Retrieve results
The results of the fitting are available as a `FitResult` structure, as returned by the `fit!` fuction.  The structure dump for the `result2` in the example above is as follows:
```julia
dump(result2)
ModelFit.FitResult
  fitter: ModelFit.Minimizer ModelFit.Minimizer()
  param: ModelFit.HashVector{ModelFit.FitParameter}
    keys: Array{Symbol}((6,))
      1: Symbol comp1__p1
      2: Symbol comp1__p2
      3: Symbol comp1__p3
      4: Symbol comp2__p1
      5: Symbol comp2__p2
      6: Symbol calib__val
    values: Array{ModelFit.FitParameter}((6,))
      1: ModelFit.FitParameter
        val: Float64 1.0034771773408524
        unc: Float64 0.0029174421450895533
      2: ModelFit.FitParameter
        val: Float64 0.0010022369285940917
        unc: Float64 1.3473265492873467e-6
      3: ModelFit.FitParameter
        val: Float64 9.996494571427457e-7
        unc: Float64 1.3148179368867407e-10
      4: ModelFit.FitParameter
        val: Float64 4.005022216534921
        unc: Float64 0.0013760967902191378
      5: ModelFit.FitParameter
        val: Float64 4.999999879104177
        unc: Float64 5.998373321161513e-8
      6: ModelFit.FitParameter
        val: Float64 1.299996982192236
        unc: Float64 4.979637254678212e-5
  ndata: Int64 399962
  dof: Int64 399956
  cost: Float64 400117.9725389123
  status: Symbol Optimal
  elapsed: Float64 1.736617397
```
From this structure the user can retrieve the parameter best fit values and uncertainties, the number of data samples, the number of degrees of freedom, the total chi-squared (`cost`) and the fitting elapsed time in seconds.

In particular, the best fit value and uncertanty for `p1` can be retrieved as follows:
```julia
println(result2.param[:comp1__p1].val)
println(result2.param[:comp1__p1].unc)
```
Note that the parameter is identified by using both the component name and parameter name, separated by a double underscore `__`.prepending the 

To plot the results use the following arrays:
- Abscissa: `dom[1]`.  For multi dimensional domains you can use `dom[2]`, `dom[3]`, etc.;
- Ordinate:
  - expression 1, i.e. `(comp1 + comp2) * comp3`: `model()` or `model2(1)`;
  - expression 2, i.e. `calib * ((comp1 + comp2) * comp3)`: `model2(2)`;
  - `comp1` component: `model2(:comp1)`;
  - `comp2` component: `model2(:comp2)`;
  - `comp3` component: `model2(:comp3)`;
  - `calib` component: `model2(:calib)`;

In the following example I will use the [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl) package, but the user can choose any other package:
```julia
using Gnuplot
@gp xr=(1000, 1010) :-
@gp :- dom[1] model2(:comp1) "w l tit 'comp1'" :-
@gp :- dom[1] model2(:comp2) "w l tit 'comp2'" :-
@gp :- dom[1] model2() "w lines tit 'Model' lw 3"
```

To evaluate the model with a different parameter value you can pass one (or more) `Pair{Symbol,Number}` to `model2()`,i.e.:
```julia
@gp :- dom[1] model2(:comp2__p1=>0.) "w lines tit 'p4=0' lw 3"
```


## The `FuncWrap` and `SimpleParam` components.

`FuncWrap` and `SimpleParam` are the only built-in components of the `ModelFit` package.  Further components are available in the `ModelFit_SimpleComponents` package.

The `FuncWrap` is simply a wrapper to a user defined function of the form `f(x, [y], [z], [further dimensions...], p1, p2, [further parameters...])`.  The number of parameters should be between 0 and 10.  If further parameters are needed the appropriate code can be generated using the `ModelFit.@code_funcwrap_npar` macro, i.e.:
```julia
ModelFit.@code_funcwrap_npar 20
```
will allow to use `FuncWrap` with 20 parameters.

The `SimpleParam` represent a scalar component in the model, whose value is given by the `val` parameter.


## Profile a component
The `ModelFit` framework provides an effective way to check a component performance by means of the `test_component` function.

For instance, you may wonder if the `FuncWrap` component used in the previous example as a wrapper to the `f` function introduces any performance penalty.  With `test_component` you may simply compare the wrapper performance with that of a simple loop, i.e.:
```julia
test_component(dom, FuncWrap(f, p1, p2, p3, p4, p5), 1000)
@time for i in 1:1000
    dummy = f(x, p1, p2, p3, p4, p5)
end
```
As you can see there is no significant difference in looping 1000 times on model evaluations or using a simple loop on `f`.



## Using the `CMPFit` minimizer

The `ModelFit` package uses the [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) minimizer by default, but it allows to use the [CMPFit](https://github.com/gcalderone/CMPFit.jl) as well.  The latter typically provides better performances with respect to the former, but since `CMPFit` is a wrapper to a C library it is not a pure-Julia solution.   The better performances are due to a different minimization strategy, not to C vs. Julia differences.

To use the `CMPFit` minimzer (after having installed the package):
```julia
result2 = fit!(model2, [data, data2], minimizer=CMPFit.Minimizer())
```


## Multidimensional domains

**IMPORTANT NOTE**: by default `the `ModelFit` package defines structure only for 1D and 2D fitting.  To handle higher dimensional cases you should generate the appropriate code with, e.g.:
```julia
ModelFit.@code_ndim 3
```

`N`-dimensional `Domain` objects comes in two flavours: linear and cartesian ones:
- Linear domain: contains a `N`-dimensional tuple of coordinates, one for each data sample.  It is analogous to `Vector{NTuple{N, Number}}`;
- Cartesian domain: contains `N` arrays of coordinates, one for each dimension. Optionally contains a 1D list of index to select a subset of all possible combinations of coordinates. It is analogous to `Vector{Vector{Number}}`, whose length of outer vector is `N`;

Linear domains are created using the `Domain` function, providing as many arguments as the number of dimensions. Each argument can either be an integer (specifying how many samples are defined along each axis), or a vector of `float`s.  **The length of all dimensions must be exactly the same.**  Examples:
```julia
# 1D
dom = Domain(5)
dom = Domain(1.:5)
dom = Domain([1,2,3,4,5.])

# 2D
dom = Domain(5, 5)
dom = Domain(1.:5, [1,2,3,4,5.])
```
Note that the length of all the above domains is 5.

Cartesian domains are created using the `CartesianDomain` function, providing as many arguments as the number of dimensions.  There is no 1-dimensional cartesian domain, hence `CartesianDomain` requires at least two arguments.  Each argument can either be an integer (specifying how many samples are defined along each axis), or a vector of `float`s.  **The length of dimensions may be different.**  Examples:
```julia
# 2D
dom = CartesianDomain(5, 6)
dom = CartesianDomain(1.:5, [1,2,3,4,5,6.])
```
The length of all the above domains is 30, i.e. it is equal to the product of the lengths of all dimensions.

Typically, the model can be evaluated over the cartesian product of all dimensions.  However, there can be specific locations of the domain for which there is not empirical data to compare with, making the model evaluation useless.  In these cases it is possible to select a subset of the cartesian domain using a 1D linear index, e.g.:
```julia
dom = CartesianDomain(1.:5, [1,2,3,4,5,6.], index=collect(0:4) .* 6 .+ 1)
```
The length of such domain is 5, equal to the length of the vector passed as `index` keyword.


A cartesian domain can always be transformed into a linear domain, while the inverse is usually not possible.  To check how a "flattened" version of a cartesian domain looks like you can use the `ModelFit.flatten` function, i.e.:
```julia
dom = CartesianDomain(1.:5, [1,2,3,4,5,6.], index=collect(0:4) .* 6 .+ 1)
lin = ModelFit.flatten(dom)
```
Actually, all models are always evaluated on "flattened", i.e. linear, domains.

To get the vector of coordinates for dimensions 1, 2, etc. of the `dom` object use the `dom[1]`, `dom[2]`, etc. syntax.  For linear domains all such vectors have the same length, for cartesian domains the lengths may differ.


## Multidimensional fitting

As an example we will fit a 2D plane.  The analytic function will accept two vectors for coordinates x and y, and 2 parameters.
```julia
f(x, y, p1, p2) = @.  p1 * x  +  p2 * y
```
This function will be called during model evaluation, where only linear domains are involved, hence we are sure that both the `x` and `y` vectors will have the same lengths.

To create a multidimensional `Measures` object we will populate a 2D array and use it as argument to the `Measures` function:
```julia
dom = CartesianDomain(30, 40)
d = fill(0., size(dom));
for i in 1:length(dom[1])
    for j in 1:length(dom[2])
        d[i,j] = f(dom[1][i], dom[2][j], 1.2, 2.4)
    end
end
data = Measures(d + randn(rng, size(d)), 1.)
```

To fit the model proceed in the usual way:
```julia
model = Model(:comp1 => FuncWrap(f, 1, 2))
prepare!(model, dom, :comp1)
result = fit!(model, data)
```


## Parameter settings
Each model parameter has a few settings that can be tweaked by the user before running the actul fit:
- `.val` (float): guess initial value;
- `.low` (float): lower limit for the value (default: `-Inf`);
- `.high` (float): upper limit for the value (default: `+Inf`);
- `.fixed` (bool): false for free parameters, true for fixed ones (default: `false`);
- `.expr` (string): a mathematical expression to bind the parameter value to other parameters (default: `""`);

**Important note**: the default minimizer ([LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) do not supports bounded parameters, while  [CMPFit](https://github.com/gcalderone/CMPFit.jl)) supports them.

Considering the previous example we can limit the interval for `p1`, and fix the value for `p2` as follows:
```julia
model = Model(:comp1 => FuncWrap(f, 1, 2))
model.comp[:comp1].p1.val  = 1   # guess initial value
model.comp[:comp1].p1.low  = 0.5 # lower limit
model.comp[:comp1].p1.high = 1.5 # upper limit
model.comp[:comp1].p2.val  = 2.4
model.comp[:comp1].p2.fixed = true
prepare!(model, dom, :comp1)
result = fit!(model, data, minimizer=CMPFit.Minimizer())
```


As another example we may constrain `p2` to always be twice the value of `p1`:
```julia
model = Model(:comp1 => FuncWrap(f, 1, 2))
model.comp[:comp1].p2.expr = "2 * comp1__p1"
model.comp[:comp1].p2.fixed = true
prepare!(model, dom, :comp1)
result = fit!(model, data)
```
Note that we had to fix the `p2` parameter otherwise the minizer will try to find a best fit for a parameter which has no influence on the final model, since its value will always be overwritten by the expression.  The only situation where the parameter is to be left free is when the expression involve the parameter itself, e.g.:
```julia
model.comp[:comp1].p2.expr = "comp1__p1 + comp1__p2"
```


## Component settings
A model component can be disabled altoghether and its evaluation fixed to a specific value  by calling `setcompvalue!`, e.g.:
```julia
setcompvalue!(model, :comp1, 10.)
```
In further evaluation the `comp1` component will always evaluate to 10.  Note thet the component parameters are completely ignored in this case, so be sure to have at least one free parameter before running the fit.

To restore the normal component behaviour simply set its value to `NaN`, i.e.:
```julia
setcompvalue!(model, :comp1, NaN)
```
