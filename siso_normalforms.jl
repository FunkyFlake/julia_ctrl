using ControlSystems, LinearAlgebra

# Let x = T * x_transformed
function similaritySS(T, T_inv, A, B, C, D) 
    A_sim = T_inv * A * T
    B_sim = T_inv * B
    C_sim = C * T
    D_sim = D
    return A_sim, B_sim, C_sim, D_sim
end

# General SISO system
A = [0  0  0.5; 
    -2 -2 -1.5;
    -6 -4 -3]
B = [0; 
     1; 
     2]
C = [3 2 -0.5]
D = 0

# Regelungsnormalform
Qs = ctrb(A,B)
Qs_inv = inv(Qs)
e_n = Qs_inv[end, :]' # last row of Qs_inv
# multiply e_n with each power of A 
T_inv = reduce(vcat, [e_n * A^i for i ∈ 0:size(A)[1]-1])
T = inv(T_inv)

A_rnf, B_rnf, C_rnf, D_rnf = similaritySS(T, T_inv, A, B, C, D)
A_rnf
B_rnf
C_rnf
D_rnf

# Beobachtungsnormalform
Qb = obsv(A,C)
Qb_inv = inv(Qb)
e_n = Qb_inv[:, end] # last column of Qb_inv
# multiply e_n with each power of A
T = reduce(hcat,[A^i * e_n for i ∈ 0:size(A)[1]-1])
T_inv = inv(T)

A_bnf, B_bnf, C_bnf, D_bnf = similaritySS(T, T_inv, A, B, C, D)
A_bnf
B_bnf
C_bnf
D_bnf


