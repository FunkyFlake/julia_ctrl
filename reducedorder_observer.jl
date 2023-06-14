using ControlSystems, LinearAlgebra
using Plots

# Observer with s states measured directly -> only order n-s observer needed
n = 3
outs = 2 
s = 2
nobs = n-s

# Order states s.t. the measured states are the first s states
A = [0  0   0.4989    0       0;
     0  0   0         4.008   0;
     0  0  -0.4022    0       0;
     0  0  -0.009793 -0.7637 -0.9658;
     0  0   0.03684  -1.498  -3]
B = [ 0      0;
      0      0;
     12.04   0;
     -1.728 -6.541;
      6.680 -6.984]
C = [1 0 0 0 0;
     0 1 0 0 0] # [I | 0] Identity measurement matrix
D = [0 0;
     0 0]
sys = ss(A, B, C, D)

# Reduced observer requires full observability of (A,C)
rank(obsv(A,C))

# Pole placement
# Desired poles for observer of unmeasured states
A12 = A[1:s,s+1:end]     # "C"
A22 = A[s+1:end,s+1:end] # "A" 
B2 = B[s+1:end,:]        # "B"

Q = [2 0 0;
     0 2 0;
     0 0 1]
R = [0.1 0;
     0 0.1]

# Substitution of A->A', B->C' and K->L' for observer design of output feedback L(y-yhat)
# For guaranteed exponential convergence rate e^-αt, α>0: all poles of observer left of -α
α = 0
A_lqr = A22' + α*I
B_lqr = A12'
# LQR requires (A_lqr,Q1) to be observable where Q = Q1'Q1
# algebraic riccati equation A'P + PA - PBR^{-1}B'P + Q = 0 
# where A refers to A_lqr and B to B_lqr
P = are(Continuous, A_lqr, B_lqr, Q, R)
# Ideal feedback gain K = R^{-1}B'P
K = R \ B_lqr'P

# Direct calculation of the feedback gain K
K = lqr(A_lqr, B_lqr, Q, R)

L = K'

# System with uncertainties
Am = A + I*0.1
Bm = B .+ 0.1
Cm = C
Dm = D

# Build system with observer [xdot xhatdot]' = [A 0; LC Am-LCm] [x xhat]' + [B; Bm] u
Aobs = [A      zeros(size(A));
        L*C          Am-L*Cm]
Bobs = [B;
        Bm]
Cobs = [C               zeros(size(C));
        zeros(size(C))              Cm]
D_obs = [D;
         Dm]

sys_obs = ss(Aobs, Bobs, Cobs, D_obs)

# Simulate
x0 = [0; 1; -1]
xhat0 = [0.5; 0.5; 0.5]
t = 0:0.01:5
u(x,t) = [sin(t);
          0.5*cos(t)]
yobs, t, xobs = lsim(sys_obs, u, t, [x0; xhat0])
# yobs, t, xobs = lsim(sys_obs, zeros(outs,length(t)), t, [x0; xhat0])
y = yobs[1:outs,:]
yhat = yobs[outs+1:end,:]
x = xobs[1:n,:]
xhat = xobs[n+1:end,:]
e = y - yhat # observer error

# Plot
# Create subplots
# set legend to the right
p1 = plot(t,[x' xhat'], lab=["x1" "x2" "x3" "xhat1" "xhat2" "xhat3"], legend=:topright, title="State and Observer states")
p2 = plot(t,[y' yhat'], lab=["y1" "y2" "yhat1" "yhat2"], title="Output and Observer output")
p3 = plot(t,e', lab=["e1" "e2"], title="Observer error")
plot(p1, p2, p3, layout=(3,1), size=(800,600))
