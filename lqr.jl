using ControlSystems, LinearAlgebra, Plots

A = [0   1;
     -2 -5]
B = [1;
     2]
C = [1 0]
D = [0]
sys = ss(A, B, C, D)

Q = [1 0;
     0 2]    
R = 1I

# For guaranteed exponential convergence rate e^-αt, α>0: all poles of A-BK left to -α
α = 2
A_lqr = A + α*I
# LQR requires (A,Q1) to be observable where Q = Q1'Q1
# algebraic riccati equation A'P + PA - PBR^{-1}B'P + Q = 0
P = are(Continuous, A_lqr, B, Q, R)
# Ideal feedback gain K = R^{-1}B'P
K = R \ B'P

# Direct calculation of the feedback gain K
K = lqr(A_lqr, B, Q, R)

# Closed loop system
A_cl = A - B*K
B_cl = B
C_cl = C
D_cl = D
sys_cl = ss(A_cl, B_cl, C_cl, D_cl)

# Simulate
t = 0:0.01:10
x0 = [1; 1]

# Open Loop
y, t, x = lsim(sys, zeros(length(t))', t, x0)
plot(t,x', lab=["x1" "x2"])
u(x, t) = -K*x
y_cl, t, x_cl, uout = lsim(sys,u,t,x0=x0)
plot!(t,x_cl', lab=["x1" "x2"])
plot!(t,uout', lab="u")