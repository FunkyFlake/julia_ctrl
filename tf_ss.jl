using ControlSystemsBase, Plots

tf(1, [1,2])

s = tf("s")

P = (s+0.5) / (s^2 + 2s + 1)
C = pid(1,2)

numvec(P)
denvec(P)
numpoly(P)

2P
inv(P)
P+C
P*C

[P; C]
[P C]

P = ss(P)
P.nx # number of states
P.nu # number of inputs
P.ny # number of outputs

P*C # ss * tf

G1 = P*C / (1 + P*C)
G2 = feedback(P*C, 1) # vorwaerts, rueckwaerts

A = [0 1; -1 -2]
B = [0; 1]
C = [1 0]
D = 0
P = ss(A, B, C, D)

P.A
(; A, B, C, D) = P

poles(P)
using LinearAlgebra
eigvals(P.A)

tzeros(P)
Pd = c2d(P, 0.1) # 0.1 sampling time, default: zero-order hold
Pd.Ts

ControlSystemsBase.timeevol(Pd) # time evolution

Ts = 0.1
A = [1 Ts; 0 1]
Pd2 = ss(A, B, C, D, Ts) # discrete statespace

isdiscrete(Pd)

ControlSystemsBase.issiso(Pd)
isstable(Pd)

C = pid(1,2)
c2d(C, 0.1)


ps = [-1, -2]
zs = [-3]
P = zpk(zs, ps, 1) # zero-pole gain

using Plots
bodeplot(P)

ω = exp10.(LinRange(-2,2,200))
bodeplot(P, ω, c=:red, plotphase=false, label="P", title="Bode Plot")