"""juliaでやる"""


using Plots
using LinearAlgebra

using StaticArrays
using ArraysOfArrays


function split_vec_of_arrays(u)
    vec.(u) |>
    x ->
    VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end

function solve_RungeKutta(dx, x₀::Vector{T}, t_span, Δt::T) where T
    """ルンゲクッタ法（4次）"""

    t = range(t_span..., step = Δt)  # 時間軸
    x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

    x[1] = x₀  # 初期値
    for i in 1:length(x)-1
        k₁ = dx(x[i])
        k₂ = dx(x[i]+k₁*Δt/2)
        k₃ = dx(x[i]+k₂*Δt/2)
        k₄ = dx(x[i]+k₃*Δt)
        x[i+1] = x[i] .+ (k₁ .+ 2k₂ .+ 2k₃ .+k₄) .* Δt/6
    end

    t, x
end

include("derived/N_is_1/eqs/julia_style/C0_0.jl")
include("derived/N_is_1/eqs/julia_style/C0_1.jl")
include("derived/N_is_1/eqs/julia_style/C0_2.jl")
include("derived/N_is_1/eqs/julia_style/C1_0.jl")
include("derived/N_is_1/eqs/julia_style/C1_1.jl")
include("derived/N_is_1/eqs/julia_style/C1_2.jl")
include("derived/N_is_1/eqs/julia_style/C2_0.jl")
include("derived/N_is_1/eqs/julia_style/C2_1.jl")
include("derived/N_is_1/eqs/julia_style/C2_2.jl")

include("derived/N_is_1/eqs/julia_style/M0_0.jl")
include("derived/N_is_1/eqs/julia_style/M0_1.jl")
include("derived/N_is_1/eqs/julia_style/M0_2.jl")
include("derived/N_is_1/eqs/julia_style/M1_0.jl")
include("derived/N_is_1/eqs/julia_style/M1_1.jl")
include("derived/N_is_1/eqs/julia_style/M1_2.jl")
include("derived/N_is_1/eqs/julia_style/M2_0.jl")
include("derived/N_is_1/eqs/julia_style/M2_1.jl")
include("derived/N_is_1/eqs/julia_style/M2_2.jl")

include("derived/N_is_1/eqs/julia_style/G0.jl")
include("derived/N_is_1/eqs/julia_style/G1.jl")
include("derived/N_is_1/eqs/julia_style/G2.jl")

include("derived/N_is_1/eqs/julia_style/Phi0.jl")
include("derived/N_is_1/eqs/julia_style/Theta0.jl")


using .C0_0
using .C0_1
using .C0_2
using .C1_0
using .C1_1
using .C1_2
using .C2_0
using .C2_1
using .C2_2

using .M0_0
using .M0_1
using .M0_2
using .M1_0
using .M1_1
using .M1_2
using .M2_0
using .M2_1
using .M2_2

using .G0
using .G1
using .G2
using .Phi0
using .Theta0


const N = 1  # セクションの数


"""慣性行列"""
function M(q_large::Vector{T}, xi_large::Vector{T}, q_dot_large::Vector{T}) where T
    Z = Matrix{T}(undef, 3*N, 3*N)

    Z[1,1] = M0_0.f(q_large, xi_large, q_dot_large)
    Z[1,2] = M0_1.f(q_large, xi_large, q_dot_large)
    Z[1,3] = M0_2.f(q_large, xi_large, q_dot_large)
    Z[2,1] = M1_0.f(q_large, xi_large, q_dot_large)
    Z[2,2] = M1_1.f(q_large, xi_large, q_dot_large)
    Z[2,3] = M1_2.f(q_large, xi_large, q_dot_large)
    Z[3,1] = M2_0.f(q_large, xi_large, q_dot_large)
    Z[3,2] = M2_1.f(q_large, xi_large, q_dot_large)
    Z[3,3] = M2_2.f(q_large, xi_large, q_dot_large)

    return Z
end


"""コリオリ＋遠心力行列"""
function C(q_large::Vector{T}, xi_large::Vector{T}, q_dot_large::Vector{T}) where T
    Z = Matrix{T}(undef, 3*N, 3*N)

    Z[1,1] = C0_0.f(q_large, xi_large, q_dot_large)
    Z[1,2] = C0_1.f(q_large, xi_large, q_dot_large)
    Z[1,3] = C0_2.f(q_large, xi_large, q_dot_large)
    Z[2,1] = C1_0.f(q_large, xi_large, q_dot_large)
    Z[2,2] = C1_1.f(q_large, xi_large, q_dot_large)
    Z[2,3] = C1_2.f(q_large, xi_large, q_dot_large)
    Z[3,1] = C2_0.f(q_large, xi_large, q_dot_large)
    Z[3,2] = C2_1.f(q_large, xi_large, q_dot_large)
    Z[3,3] = C2_2.f(q_large, xi_large, q_dot_large)

    return Z
end


"""重力行列"""
function G(q_large::Vector{T}, xi_large::Vector{T}, q_dot_large::Vector{T}) where T
    Z = Vector{T}(undef, 3*N)

    Z[1] = G0.f(q_large, xi_large, q_dot_large)
    Z[2] = G1.f(q_large, xi_large, q_dot_large)
    Z[3] = G2.f(q_large, xi_large, q_dot_large)

    return Z
end


"""状態方程式"""
function X_dot(X::Vector{T}) where T
    q = X[1:3]
    q_dot = X[4:6]
    xi = [1.0, 1.0, 1.0]

    τ = [0.0, 0.0, 0.0]

    inv_M = inv(M(q, xi, q_dot))

    x1_dot = q_dot
    x2_dot = -inv_M * C(q, xi, q_dot) * q_dot .+
    inv_M * (τ .- G(q, xi, q_dot))

    println([x1_dot; x2_dot])
    return [x1_dot; x2_dot]

end


function test()

    TIME_SPAN = 10.0

    q_large = [0.01, 0.02, 0.0]
    xi_large = [1.0, 1.0, 1.0]
    q_dot_large = [0.0, 0.0, 0.0]

    t, X = solve_RungeKutta(
        X_dot,
        [q_large; q_dot_large],
        (0.0, TIME_SPAN),
        0.001
    )

    l1, l2, l3, l1_dot, l2_dot, l3_dot = split_vec_of_arrays(X)
    fig = plot(xlim=(0, TIME_SPAN))
    plot!(fig, t, l1, label="l1")
    plot!(fig, t, l2, label="l2")
    plot!(fig, t, l3, label="l3")
    

    fig2 = plot(xlim=(0, TIME_SPAN))
    plot!(fig2, t, l1_dot, label="l1_dot")
    plot!(fig2, t, l2_dot, label="l2_dot")
    plot!(fig2, t, l3_dot, label="l3_dot")

    fig_I = plot(fig, fig2, layout=(1, 2), size=(800, 400))
    savefig(fig_I, "julia.png")



    # アニメ制作



end


@time test()