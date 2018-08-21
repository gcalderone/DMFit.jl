if VERSION < v"0.7.0"
    using Base.Test
else
    using Test
    using Random
end
rng = MersenneTwister(0);
using DMFit

p1 = 1
p2 = 1.e-3
p3 = 1.e-6
p4 = 4
p5 = 5

f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)
x = 1:0.05:10000
y = f(x, p1, p2, p3, p4, p5);
noise = randn(rng, length(x));


dom = Domain(x)
data = Measures(y + noise, 1.)

model1 = Model(:comp1 => FuncWrap(f, p1, p2, p3, p4, p5))
prepare!(model1, dom, :comp1)

model1[:comp1].p1.val = 1
model1[:comp1].p2.val = 1.e-3

result1 = fit!(model1, data)


test_component(dom, FuncWrap(f, p1, p2, p3, p4, p5), 1000)
@time for i in 1:1000
    dummy = f(x, p1, p2, p3, p4, p5)
end


f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @.  p4 * sin.(p5 * x)
f3(x) = cos.(x)

model2 = Model(:comp1 => FuncWrap(f1, p1, p2, p3),
               :comp2 => FuncWrap(f2, p4, p5),
               :comp3 => FuncWrap(f3))
prepare!(model2, dom, :((comp1 + comp2) * comp3))
result2 = fit!(model2, data)




y2 = f(x, p1, p2, p3, p4, p5);
noise = randn(rng, length(x));

data2 = Measures(1.3 .* (y + noise), 1.3)

push!(model2, :calib, SimpleParam(1))
prepare!(model2, dom, :(calib * ((comp1 + comp2) * comp3)))
result2 = fit!(model2, [data, data2])




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


model = Model(:comp1 => FuncWrap(f, 1, 2))
model.comp[:comp1].p1.val  = 1   # guess initial value
model.comp[:comp1].p1.low  = 0.5 # lower limit
model.comp[:comp1].p1.high = 1.5 # upper limit
model.comp[:comp1].p2.val  = 2.4
model.comp[:comp1].p2.fixed = true
prepare!(model, dom, :comp1)
result = fit!(model, data, minimizer=CMPFit.Minimizer())


model = Model(:comp1 => FuncWrap(f, 1, 2))
model.comp[:comp1].p2.expr = "2 * comp1__p1"
model.comp[:comp1].p2.fixed = true
prepare!(model, dom, :comp1)
result = fit!(model, data)
