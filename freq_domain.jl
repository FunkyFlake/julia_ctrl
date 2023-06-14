using ControlSystemsBase, Plots

P = DemoSystems.double_mass_model()

pzmap(P,strokewidth=5,markersize=10)

dampreport(P)

using ControlSystems # includes more e.g. root locus

rlocusplot(P, K=10000) # K: max gain

bodeplot(P)

w = exp10.(LinRange(-2,2,200))

Pw = freqresp(P, w) # array [#inputs, #outputs, #freqs] = [1,1,200] (useful for MIMO)

Pw = freqrespv(P, w) # output as vector

# Mimo example with random ss system 
Pw = freqresp(ssrand(2,3,4), w) # 2 outputs, 3 inputs, 4 states

mag, phase = bode(P, w)

bodeplot(P, w, title="Double Mass Model")

re, im = nyquist(P, w)
nyquistplot(P, w, Ms_circles=[2], ylims=(-4,1), xlims=(-3,1))

# Interactive plot: execute plotly() to switch plotting backend


# MIMO example
P = ssrand(2, 3, 4, proper=true)
# 2 outputs, 3 inputs, 4 states -> 2x3 transfer function matrix
# we have min(inputs, outputs) singular values, these are plotted
 
bodeplot(P, w)

# Singular value plot
sigmaplot(P, w)

dcgain(P) # per input-output pair

## Gang of four
P = DemoSystems.double_mass_model()
C = pid(10, 10, 1, Tf=0.001) # pid with 2nd order lowpass
S, PS, CS, T = gangoffour(P, C)
# S: sensitivity function 1/(1+PC)
# PS: load disturbance to measurement signal P/(1+PC)
# CS: measurement noise to control signal C/(1+PC)̧ 
# T: complementary sensitivity function PC/(1+PC)
plot(step(T, 0:0.01:10))

gangoffourplot(P,C, w, label="")

# Relative gain array
s = tf("s")
P = [1/(s+1) 1/(s+2); 1/(s+3) 1/(s+4)]

relative_gain_array(P, 0)

abs.(relative_gain_array(P, 0))

rgaplot(P, w)

# Frequency-response data
using ControlSystemIdentification
ω = 1
ξ = 0.1
G = tf(ω^2, [1, 2ξ*ω, ω^2])

frd = FRD(w, G)
plot(frd)
