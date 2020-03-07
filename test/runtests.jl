using Random

f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)

# Actual model parameters:
params = [1, 1.e-3, 1.e-6, 4, 5]

# Domain for model evaluation
x = 1.:50:10000

# Evaluated model
y = f(x, params...);

# Random noise
rng = MersenneTwister(0);
noise = randn(rng, length(x));

using DataFitting, DataFitting.Components
data = Measures(y + noise, 1.)

model = Model(:comp1 => FuncWrap(f, params...))
add_dom!(model, x)
addexpr!(model, 1, :comp1)

model.comp1.p[1].val = 1
model.comp1.p[2].val = 1.e-3

setparamvalues!(model, params)
result = fit!(model, data)

using CMPFit
DataFitting.@enable_CMPFit
setparamvalues!(model, params)
result = fit!(model, data, minimizer=cmpfit())




addcomp!(model, :ss => DataFitting.Smooth(2))
addexpr!(model, 1, :(ss(comp1)))



f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @. p4 * sin(p5 * x)
f3(x) = cos.(x)

model = Model(:comp1 => FuncWrap(f1, params[1], params[2], params[3]),
              :comp2 => FuncWrap(f2, params[4], params[5]),
              :comp3 => FuncWrap(f3))
add_dom!(model, x)
addexpr!(model, :((comp1 .+ comp2) .* comp3))
result = fit!(model, data)



p = probe(model.comp1.p[1], data)
@gp p[:,1] p[:,2] "w l notit"

p = probe([model.comp1.p[1], model.comp2.p[2]], data)
@gsp p[:,1] p[:,2] p[:,3] "w p lc palette"

@gsp "set contour" :-
@gsp :- "set dgrid3d 11,11" :-
@gsp :- "set cntrparam lev incremental 0, 1, 10" :-
@gsp :- p[:,1] p[:,2] p[:,3] "w l notit lc palette"


p = DataFitting.probe([model.comp1.p[1], model.comp2.p[2], model.comp1.p[2]], data, nstep=[11,5,5])
c = p[:,4]
c ./= maximum(c)
c = c .^ 0.5

color = Int.(
    2^24 .* round.(255 .* c                            ) .+ 
    2^16 .* round.(255 .* abs.(sin.(pi   .* c .+ pi/2))) .+
    2^08 .* round.(255 .* abs.(sin.(pi   .* c .- pi/4))) .+ 
    2^00 .* round.(255 .* abs.(sin.(pi/4 .* c .+ pi/4)))  )
@gsp p[:,1] p[:,2] p[:,3] 1.5 .*(c.+0.1) color " w p notit pt 7 ps var lc rgb var"




noise = randn(rng, length(x));
data2 = Measures(1.3 * (y + noise), 1.3)

addcomp!(model, :calib=>ScalarParam(1))
add_dom!(model, x)
addexpr!(model, 2, :(calib .* ((comp1 .+ comp2) .* comp3)))
result = fit!(model, [data, data2])



println(result.comp1.p[1].val)
println(result.comp1.p[1].unc)


test_component(FuncWrap(f, params...), x; iter=1000)
@time for i in 1:1000
    dummy = f(x, params...)
end

@eval DataFitting @code_ndim 4

# 1D
dom = Domain(5)
dom = Domain(1.:5)
dom = Domain([1,2,3,4,5.])

# 2D
dom = Domain(5, 5)
dom = Domain(1.:5, [1,2,3,4,5.])

# 2D
dom = CartesianDomain(5, 6)
dom = CartesianDomain(1.:5, [1,2,3,4,5,6.])

dom = CartesianDomain(1.:5, [1,2,3,4,5,6.], index=collect(0:4) .* 6 .+ 1)

dom = CartesianDomain(1.:5, [1,2,3,4,5,6.], index=collect(0:4) .* 6 .+ 1)
lin = DataFitting.flatten(dom)






f(x, y, p1, p2) = @.  p1 * x  +  p2 * y

dom = CartesianDomain(30, 40)
d = fill(0., size(dom));
for i in 1:length(dom[1])
    for j in 1:length(dom[2])
        d[i,j] = f(dom[1][i], dom[2][j], 1.2, 2.4)
    end
end
data = Measures(d + randn(rng, size(d)), 1.)

model = Model(:comp1 => FuncWrap(f, 1, 2))
add_dom!(model, dom)
addexpr!(model, :comp1)
result = fit!(model, data)


model.comp1.p[1].val  = 1   # guess initial value
model.comp1.p[1].low  = 0.5 # lower limit
model.comp1.p[1].high = 1.5 # upper limit
model.comp1.p[2].val  = 2.4
model.comp1.p[2].fixed = true
result = fit!(model, data)



model.comp1.p[1].low  = -Inf
model.comp1.p[1].high = +Inf


model.comp1.p[2].expr = "2 * comp1_p1"
model.comp1.p[2].fixed = true
result = fit!(model, data)

model.comp1.p[2].expr = "comp1_p1 + comp1_p2"
model.comp1.p[2].fixed = false
result = fit!(model, data)



rng = MersenneTwister(0);


dom = Domain(1:0.001:10)
test_component(DataFitting.Components.OffsetSlope(1.1, 2.2, 3.3), dom, iter=1000);

dom = CartesianDomain(100, 100)
test_component(DataFitting.Components.OffsetSlope(1.1, 2.2, 3.3, 4.4, 5.5), dom, iter=1000);

test_component(DataFitting.Components.Polynomial(1.1, 2.2, 3.3), 1:0.001:10, iter=1000);
test_component(DataFitting.Components.Gaussian(1.1 , 4.4, 0.51), 1:0.001:10, iter=1000);

dom = CartesianDomain(100, 100)
test_component(DataFitting.Components.Gaussian(100, 30, 70, 5), dom, iter=1000);
test_component(DataFitting.Components.Gaussian(100, 30, 70, 5, 0.5, 10), dom, iter=1000);


dom = Domain(1:0.001:10)
test_component(DataFitting.Components.Lorentzian(1.1 , 4.4, 0.51), dom, iter=1000);

dom = CartesianDomain(100, 100)
test_component(DataFitting.Components.Lorentzian(100, 30, 70, 5, 0.5), dom, iter=1000);


x = Domain(1:0.05:10)
model = Model(
    :offset => ScalarParam(4),
    :line1  => DataFitting.Components.Gaussian(1.1 , 4.4, 0.51),
    :line2  => DataFitting.Components.Gaussian(0.52, 5.5, 1.2 ))
add_dom!(model, x)
add_expr!(model, :(offset + line1 + line2))

noise = maximum(model()) * 0.01
data = Measures(model() + noise * randn(rng, length(x)), noise);
ret1 = fit!(model, data)


model = Model()
add_comp!(model, :background => DataFitting.Components.OffsetSlope(0, 0, 0., 2., 3.))
add_comp!(model, :psf        => DataFitting.Components.Gaussian(100., 0., 0., 1, 0.3, 15))
dom = CartesianDomain(-5:0.1:5, -4:0.1:4)
add_dom!(model, dom)
add_expr!(model, :(background + psf))

noise = maximum(model()) * 0.1
data = Measures(model() .+ 4 .+ noise .* randn(length(dom)), noise);
ret1 = fit!(model, data)
