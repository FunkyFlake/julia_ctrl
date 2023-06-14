using ControlSystemsBase

A = [-1.0 1.0; 0.0 -1e-6]
B = [1.0; 0.0;;]
C = [1.0 0.0]
D = [0.0;;]

P = ss(A,B,C,D)

@show poles(P)
@show tzeros(P)

# Controllability and gramian
using LinearAlgebra
Qc = ctrb(A,B)
rank(Qc) == P.nx 

Qo = obsv(A,C)
rank(Qo)
# Checking the rank is numerically hard -> check gramian instead
 
# Ctrb gramian = ∫0∞ exp(A*t)*B*B'*exp(A'*t) dt
Wc = gram(P, :c)

# Obsv gramian = ∫0∞ exp(A'*t)*C'*C*exp(A*t) dt
Wo = gram(P, :o)

@show svdvals(Wc) # singular values
@show svdvals(Wo) 

using RobustAndOptimalControl
hsvd(P) # hankel singular values (svd of √Wc*Wo), measure of energy per state

Pr, G, T = balreal(P) # balanced realization (minimal realization)
# Pr is the reduced order system, G is gramian of reduced system, T is similarity transform matrix
# obsv and ctrb gramians are now equal
gram(Pr, :c) == gram(Pr, :o)

# equivalent: 
minreal(P)

## System with delay
W = DemoSystems.woodberry()

using Plots
plot(
    bodeplot(W, label=""),
    nyquistplot(W, xlims=(-3, 12), ylims=(-7, 3), titlefontsize=9),
    labelfontsize=9,
    size=(1500, 700))

relative_gain_array(W, 0)
rgaplot(W)

# Algebraic ricatti equation (are)
Q = I(2)
R = I(1)
X = are(P, Q, R) # no care/dare since continuous/discrete is encoded in system Plots

L = lqr(P, Q, R)

X = are(Continuous, A, B, Q, R) # continuous/discrete necessary whenn passing raw matrices

Σ = lyap(P, Q) # solves lyapunov equation for P (continuous inferred from P), x' P x is then a lyapunov equation

## Margins
P = 100DemoSystems.double_mass_model() * tf(1, [0.002, 1])
w = exp10.(LinRange(-2,2,100))
marginplot(P, w)

delaymargin(P)

# gangoffour
C = pid(1, 1, 1, Tf=0.002)
gangoffourplot(P, C, lab="", size=(1500,700))

Gcl = extended_gangoffour(P, C)
bodeplot(Gcl, lab=["S" "CS" "PS" "T"], plotphase=false, size=(1500,700))

# find package:
parentmodule(extended_gangoffour)

S = sensitivity(P, C)
Ms, ωMs = hinfnorm2(S) # robust but slower version of hinfnorm from ControlSystems

f1 = bodeplot(Gcl, lab=["S" "CS" "PS" "T"], plotphase=false, size=(1500,700))
hline!([Ms], l=(:black, :dash), linewidth=20)
scatter!([ωMs], [Ms])

w = exp10.(LinRange(-2,2,2000))
f2 = nyquistplot(P * C, w, Ms_circles = Ms)
plot(f1,f2, size=(1800,700))

dampreport(P)

# Diskmargin
parentmodule(diskmargin)
dm = diskmargin(P)
plot(dm)


## Passivity
passivityplot(P)
passivityplot(tf(1, [1, 1]))

ncfmargin(P, C)