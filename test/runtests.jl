using Random

f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)

# Actual model parameters:
params = [1, 1.e-3, 1.e-6, 4, 5]

# Domain for model evaluation
x = 1.:0.05:10000

# Evaluated model
y = f(x, params...);

# Random noise
rng = MersenneTwister(0);
noise = randn(rng, length(x));

using DataFitting
data = Measures(y + noise, 1.)

model = Model(:comp1 => FuncWrap(f, params...))
add_dom!(model, x)
addexpr!(model, 1, :comp1)

model.comp1.p[1].val = 1
model.comp1.p[2].val = 1.e-3

result = fit!(model, data)


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


p = DataFitting.probe([data], nstep=(10, 30),
                      ( model.comp1.p[1],
                        result.comp1.p[1].unc*1.3),
                      ( model.comp2.p[2],
                        result.comp2.p[2].unc*3))
@gsp p[:,1] p[:,2] p[:,3].-=minimum(p[:,3]) "w p lc palette"

@gsp "set contour" :-
@gsp :- "set dgrid3d 30,10" :-
@gsp :- "set cntrparam lev incremental 0, 1, 10" :-
@gsp :- p[:,1] p[:,2] p[:,3].-=minimum(p[:,3]) "w l lc palette"


p = DataFitting.probe([data],
                      ( model.comp1.p[1],
                       result.comp1.p[1].unc*3),
                      ( model.comp2.p[2],
                       result.comp2.p[2].unc*3),  
                      ( model.comp1.p[2], 
                       result.comp1.p[2].unc*3))
c = p[:,4].-minimum(p[:,4])



@gsp @sprintf("set cbrange [%f:%f]", extrema(c)...) :-
@gsp :- xr=extrema(p[:,1]) yr=extrema(p[:,2]) zr=extrema(p[:,3]) :-
last = 0.
for i in 0.:1.25:maximum(c)
    j = findall(last .< c .<= i)
    global last = i
    if length(j) > 1
        @gsp :- p[j,1] p[j,2] p[j,3] c[j] " u 1:2:3:(1./\$4):4 w p notit pt 1 ps var lc palette"
        sleep(0.05)
    end
end








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



