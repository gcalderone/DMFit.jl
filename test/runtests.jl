using Test
using Random


f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)

# Actual model parameters:
params = [1, 1.e-3, 1.e-6, 4, 5]

# Domain for model evaluation
x = 1:0.05:10000

# Evaluated model
y = f(x, params...);

# Random noise
rng = MersenneTwister(0);
noise = randn(rng, length(x));

using DataFitting
dom = Domain(x)
data = Measures(y + noise, 1.)

model1 = Model(:comp1 => FuncWrap(f, params...))
prepare!(model1, dom, :comp1)

model1.comp1.p[1].val = 1
model1.comp1.p[2].val = 1.e-3

result1 = fit(model1, data)


f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @. p4 * sin(p5 * x)
f3(x) = cos.(x)

model2 = Model(:comp1 => FuncWrap(f1, params[1], params[2], params[3]),
               :comp2 => FuncWrap(f2, params[4], params[5]),
               :comp3 => FuncWrap(f3))
prepare!(model2, dom, :((comp1 + comp2) * comp3))
result2 = fit(model2, data)



noise = randn(rng, length(x));
data2 = Measures(1.3 * (y + noise), 1.3)

push!(model2, :calib=>SimpleParam(1))
prepare!(model2, dom, :(calib * ((comp1 + comp2) * comp3)))
result2 = fit(model2, [data, data2])




resetcounters!(model2)


dump(result2)

println(result2.bestfit.comp1.p1.val)
println(result2.bestfit.comp1.p1.unc)


test_component(dom, FuncWrap(f, params...), 1000)
@time for i in 1:1000
    dummy = f(x, params...)
end

@eval DataFitting @code_ndim 3

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
prepare!(model, flatten(dom), :comp1)
result = fit(model, data)


model.comp1.p[1].val  = 1   # guess initial value
model.comp1.p[1].low  = 0.5 # lower limit
model.comp1.p[1].high = 1.5 # upper limit
model.comp1.p[2].val  = 2.4
model.comp1.p[2].fixed = true
result = fit(model, data)



model.comp1.p[1].low  = -Inf
model.comp1.p[1].high = +Inf


model.comp1.p[2].expr = "2 * comp1__p1"
model.comp1.p[2].fixed = true
prepare!(model)
result = fit(model, data)

model.comp[:comp1].p[2].expr = "comp1__p1 + comp1__p2"
model.comp[:comp1].p[2].fixed = false
prepare!(model)
result = fit(model, data)



