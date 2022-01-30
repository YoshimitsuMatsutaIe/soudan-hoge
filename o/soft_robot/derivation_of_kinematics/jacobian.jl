using Symbolics



include("./kinematics.jl")
using .Kinematics


N = 2
@variables Q[1:3N]
@variables Ξ[1:N]

println(Q[1])

Ps, Rs = calc_local(N, Q, Ξ);

Φs, Θs = calc_global(N, Ps, Rs);


