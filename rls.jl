using ControlSystems, LinearAlgebra
using DifferentialEquations, Plots

# Regelungstechnik 2 von Schulz, Kapitel 16.3.1 zu Recursive Least Squares

# Plant Model
T = 0.1;
J = 0.5;
F = tf([40], [J, 5, 20])

Fz = c2d(F, T, :zoh)
num = numvec(Fz)[1]           # [b1, b2]
den = denvec(Fz)[1][2:end]    # [a1 a2]
wtrue = cat(den, num; dims=1) # [a1 a2 b1 b2] ground truth
n = size(wtrue)[1]

##############################################################################
# Recursive Least Squares
mutable struct RLS
    w::Vector{Float64} # n x 1
    P::Matrix{Float64} # n x n
    α::Float64         # forgetting factor
    x::Vector{Float64} # dim(x) + dim(y) = n
    y::Vector{Float64} # depending on model
end

# ynew: output measurement of model
# xnew: input into model to prepare for next step
# finds parameters w and covariance P predicting current ynew from last ys and xs
function update!(rls::RLS, ynew::Float64, xnew::Float64)
    h = cat(rls.y, rls.x, dims=1) # build regressor
    
    rls.P /= rls.α

    rls.w = rls.w - rls.P*h*(h'*rls.w - ynew)/(1 + h'*rls.P*h) 
    rls.P = rls.P - rls.P * h * h' *rls.P' / (1 + h' * rls.P * h)
    
    # pred = rls.w' * cat(rls.y, rls.x, dims=1) # predicted output
    # print("error: $(pred-ynew)\n")

    rls.x = cat(xnew, rls.x[1:end-1]; dims=1) # shift in current x for next step
    rls.y = cat(-ynew, rls.y[1:end-1]; dims=1) # shift in new y for next step
    
    return rls.w, rls.P
end

# updates the RLS estimator with an adaptive forgetting factor that keeps the 
# trace of P constant to avoid numerical blow-up or convergence to zero
function updateTrace!(rls::RLS, ynew::Float64, xnew::Float64)
    h = cat(rls.y, rls.x, dims=1) # build regressor
    
    rls.α = tr(rls.P) / size(P)[1] 

    rls.P /= rls.α 

    rls.w = rls.w - rls.P*h*(h'*rls.w - ynew)/(1 + h'*rls.P*h) 
    rls.P = rls.P - rls.P * h * h' *rls.P' / (1 + h' * rls.P * h)
    
    # pred = rls.w' * cat(rls.y, rls.x, dims=1) # predicted output
    # print("error: $(pred-ynew)\n")

    rls.x = cat(xnew, rls.x[1:end-1]; dims=1) # shift in current x for next step
    rls.y = cat(-ynew, rls.y[1:end-1]; dims=1) # shift in new y for next step
    
    return rls.w, rls.P
end
##############################################################################

# Discrete model: y[k] = a1*y[k-1] + a2*y[k-2] + b1*x[k-1] + b2*x[k-2]
# Initialization step
w = zeros(n)                  # initial parameter estimate
P = Matrix{Float64}(I, n, n)  # initial covariance estimate
α = 0.95       # forgetting factor
x = zeros(2)   # last inputs  [x[k-1] x[k-2]] 
y = zeros(2)   # last outputs [y[k-1] y[k-2]]
rls = RLS(w, P, α, x, y)

##############################################################################
# Simulation of RLS with a plant that changes parameters
Fss = ss(F) # Plant to statespace
A, B, C, D = Fss.A, Fss.B, Fss.C, Fss.D

# Logging variables
u_in = []
w_est = []
P_est = []
y_pred = []

function f!(dx, x, p, t) 
    u = p[1]
    A,B = p[3:4]

    dx .= A*x + B*u   
    
    # Logging
    push!(u_in, u)
end

# Parameter change callback
condition2(x, t, integrator) = (t ∈ [20, 40])
function affect2!(integrator)
    if integrator.t == 20.0
        J = 2
    elseif integrator.t == 40.0
        J = 5
    end
    
    F = tf([40], [J, 5, 20])
    sys = ss(F)
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    integrator.p[3:end] = [A,B,C,D]
end
cb2 = DiscreteCallback(condition2, affect2!, save_positions=(true, true))

# Controller Callback
tspan = (0.0, 60.0)
Tc = 0.1    # controller period
tctrl = tspan[1]:Tc:tspan[2]
condition(x, t, integrator) = (t ∈ tctrl)
function affect!(integrator)
    # Calculate control signal
    if(integrator.t % 5 >= 0) && (integrator.t % 5 < 2.5)
        xnew = 1.0
    elseif (integrator.t % 5 >= 2.5) && (integrator.t % 5 < 5.0)
        xnew = -1.0
    end

    integrator.p[1] = xnew

    # Update RLS estimator
    rls = integrator.p[2]
    ynew = C*integrator.u + D*xnew              # measure output
    w,P = update!(rls, ynew[1], xnew)           # constant forgetting factor -> unstable P
    # w,P = updateTrace!(rls, ynew[1], xnew)    # not sure yet if this works correctly 
    # Logging
    push!(w_est, w)
    push!(P_est, P)
end
cb = DiscreteCallback(condition, affect!, save_positions=(true, true))

cbs = CallbackSet(cb, cb2)

x0 = zeros(size(den))
p = [1.0, rls, A,B,C,D]
prob = ODEProblem(f!, x0, tspan, p, callback=cbs, tstops=tctrl)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Evaluation
plotly()
# gr()
# System plots
plot(sol, idxs=(0, 1, 2), label=["t" "x1" "x2"], xaxis="Time [s]", yaxis="x1", zaxis="x2", size=(1500,800))
plot(sol, label=["x1" "x2"], xaxis="Time [s]", yaxis="x")
plot((C*sol(tctrl))', label="Output")
# Estimation plots
w_est = reduce(hcat,w_est)
# extract 4 vectors containing the diagonal elements of P_est
σ_a1 =[(P_est[i][1,1] >= 0 ? sqrt(P_est[i][1,1]) : 0) for i in 1:length(P_est)];
σ_a2 =[(P_est[i][2,2] >= 0 ? sqrt(P_est[i][2,2]) : 0) for i in 1:length(P_est)];
σ_b1 =[(P_est[i][3,3] >= 0 ? sqrt(P_est[i][3,3]) : 0) for i in 1:length(P_est)];
σ_b2 =[(P_est[i][4,4] >= 0 ? sqrt(P_est[i][4,4]) : 0) for i in 1:length(P_est)];

plot(tctrl[1:end-1], w_est',figsize=(800, 600), label=["a1" "a2" "b1" "b2"], xaxis="Time [s]", yaxis="Parameters", size=(1500,800))
# fill area around w_est with std deviation σ
plot!(tctrl[1:end-1], w_est', ribbon=(σ_a1, σ_a2, σ_b1, σ_b2), fillalpha=0.2, label=["σ_a1" "σ_a2" "σ_b1" "σ_b2"])
# hline!(wtrue)