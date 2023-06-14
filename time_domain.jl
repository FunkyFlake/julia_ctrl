using ControlSystemsBase, Plots

P = DemoSystems.double_mass_model(outputs=2)
P = c2d(P, 0.001) # zoh, 1ms 

res = step(P, 5) # 5s step response
plot(res)

res = impulse(P, 0:0.001:5) # step size must match sampling time
plot(res)

res.y # (num. outputs x response) matrix
res.x # state vector
plot(x', layout=4) # need to transpose to plot 4 signals with 5001 data points, otherwise 5001 signals
plot(x', layout=(4,1))
res.t
res.sys

y, t, x = impulse(P, 5)

u = randn(1, 100) # (inputs, values) matrix
res = lsim(P,u)
plot(res, plotu=true)

P = c2d(DemoSystems.double_mass_model(outputs=1:4, load_input=true), 0.001)
u = [
    randn(1,100)    # first input just random numbers
    (1:100)' .> 50  # second input, step from 0 to 1 after 50 samples
]

res = lsim(P, u)
plot(res, plotu=true, size=(600,1000), margin=5Plots.mm)

using ControlSystems # main package needed for continuous simulation
using LinearAlgebra
P = DemoSystems.double_mass_model(outputs=2)
res = step(P,5)
plot(res)

# now with continuous input u
Ts = 0.1
A = [1 Ts; 0 1]
B = [0; 1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I
R = I
L = lqr(Discrete, A,B,Q,R)

u = (x,t) -> -L*x .+ 2.5 * (t >= 2.5) # control law with disturbance after 2.5s
u([1; 2], 0) # u for two states at time 0
u([1; 2], 3) # disturbance active

t = 0:Ts:5
x0 = [1, 0]
res = lsim(sys,u,t,x0=x0)

plot(res, plotu=true)


# We can look for the implementation of library functions by typing @edit step(P,5) for example

# print function that led to a result for debugging etc. with @show macro
@show poles(P)