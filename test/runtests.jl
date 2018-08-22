if VERSION < v"0.7.0"
    using Base.Test
else
    using Test
    using Random
end


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

using DMFit
dom = Domain(x)
data = Measures(y + noise, 1.)

model1 = Model(:comp1 => FuncWrap(f, params...))
prepare!(model1, dom, :comp1)

model1[:comp1].p1.val = 1
model1[:comp1].p2.val = 1.e-3

result1 = fit!(model1, data)


f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @. p4 * sin(p5 * x)
f3(x) = cos.(x)

model2 = Model(:comp1 => FuncWrap(f1, params[1], params[2], params[3]),
               :comp2 => FuncWrap(f2, params[4], params[5]),
               :comp3 => FuncWrap(f3))
prepare!(model2, dom, :((comp1 + comp2) * comp3))
result2 = fit!(model2, data)

noise = randn(rng, length(x));
data2 = Measures(1.3 * (y + noise), 1.3)

push!(model2, :calib, SimpleParam(1))

prepare!(model2, dom, :(calib * ((comp1 + comp2) * comp3)))

result2 = fit!(model2, [data, data2])


dump(result2)

println(result2.param[:comp1__p1].val)
println(result2.param[:comp1__p1].unc)

DMFit.@code_funcwrap_npar 20

test_component(dom, FuncWrap(f, params...), 1000)
@time for i in 1:1000
    dummy = f(x, params...)
end

using CMPFit
result2 = fit!(model2, [data, data2], minimizer=CMPFit.Minimizer())

DMFit.@code_ndim 3

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
lin = DMFit.flatten(dom)




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
prepare!(model, dom, :comp1)
result = fit!(model, data)


model.comp[:comp1].p1.val  = 1   # guess initial value
model.comp[:comp1].p1.low  = 0.5 # lower limit
model.comp[:comp1].p1.high = 1.5 # upper limit
model.comp[:comp1].p2.val  = 2.4
model.comp[:comp1].p2.fixed = true
result = fit!(model, data, minimizer=CMPFit.Minimizer())



model.comp[:comp1].p1.low  = -Inf
model.comp[:comp1].p1.high = +Inf


model.comp[:comp1].p2.expr = "2 * comp1__p1"
model.comp[:comp1].p2.fixed = true
prepare!(model)
result = fit!(model, data)

model.comp[:comp1].p2.expr = "comp1__p1 + comp1__p2"
model.comp[:comp1].p2.fixed = false
prepare!(model)
result = fit!(model, data)


setcompvalue!(model, :comp1, 10)
setcompvalue!(model, :comp1, NaN)


