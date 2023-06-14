using LinearAlgebra, ControlSystems

# Reference Model
Aref = [ 0  1 0;
        -1 -2 0;
         0  0 0]
Cref = [-1 0 1;
        -2 0 2]
# Can be calculated before runtime
function wref(t)
     w1 = 1-exp(-t); 
     w2 = 2(1-exp(-t));
     return [w1; w2];
end

# System Model
A = [ 0  1   2;
     -2 -3   0;
     -2  0  -1]
B = [2 1;
     0 2;
     1 3]
C = [1  2 -1;
     2 -1  3]

n = size(A,1)  # states
o = size(C,1)  # outputs
i = size(B,2)  # inputs

# Servo Compensator
# Build N with Aref block matrix 3 times on the Diagonal
znn = zeros(size(Aref))
N = [Aref znn;
     znn  Aref] # Diagonal block matrices with same eigs as Aref
m = size(N,1)

Mi = [0 1 1]'; 
zM = zeros(size(Mi)); 
M = [Mi zM;
     zM Mi]


# Augmented state vector [x; servo]
znm = zeros(n,m)
Aaug = [ A   znm;  # normal states
        -M*C   N ] # error states
zmi = zeros(m,i)
Baug = [B;    # normal B
        zmi]  # no influence of B on e

# LQR
# Q: first 3 rows, are state 0 
Q = zeros(n+m,n+m); 
# only error is weighted
Q[4:6,4:6] .= 10;
Q[7:9,7:9] .= 1;
Q

R = [1 0;
     0 1]

K = lqr(Aaug, Baug, Q, R)
Kx = K[:,1:n]
Kservo = K[:,n+1:end]

# Simulation
using DifferentialEquations
xref0 = [1 0 1]';
x0 = zeros(n,1);
servo0 = zeros(m,1);

# Ode problem
function f!(dx, x, p, t)
    y = C*x[1:n]    # calculate output
    e = wref(t) - y # calculate error
    u = -Kx*x[1:n] -Kservo*x[n+1:end] # calculate control input
   
    dx[1:n] = A*x[1:n] + B*u            # normal states
    dx[n+1:end] = N*x[n+1:end] + M*e    # error/compensator states
    dx = x
end


tspan = (0.0, 10.0)
u0 = [x0; servo0]
prob = ODEProblem(f!, u0, tspan, wref)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8);

x = sol[1:n,1,:]
e = sol[n+1:end,1,:]
y = C*x
w = wref.(sol.t)
w = hcat(w...)

# Evaluation
using Plots
# Create Subplots
# Plot state trajectorys of x1 to x3
p1 = plot(sol.t, x', label=["x1" "x2" "x3"], xlabel="t", ylabel="x", title="State trajectorys", legend=:bottomright);
# Plot outputs and reference
p2 = plot(sol.t, y', label=["y1" "y2"], xlabel="t", ylabel="y", title="Outputs");
plot!(p2, sol.t, w', label=["w1" "w2"], legend=:bottomright);
# Plot error
p3 = plot(sol.t, e', label=["e1" "e2"], xlabel="t", ylabel="e", title="Compensator States", legend=:bottomright);

plot(p1, p2, p3, layout=(3,1), size=(800,600))