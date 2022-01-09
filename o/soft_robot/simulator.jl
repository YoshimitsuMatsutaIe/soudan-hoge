"""シミュレーション"""


using LinearAlgebra
using Plots
using ProgressBars
using Parameters

include("kinematics.jl")
include("dynamics.jl")
include("controller.jl")

using .Kinematics: Phi0, Arm
using .Dynamics: H_dot, calc_q_dot_dot
using .Controller


function split_vec_of_arrays(u)
    vec.(u) |>
    x ->
    VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end


"""ルンゲクッタ法（4次）"""
function solve_RungeKutta(dx, x₀::Vector{T}, t_span, Δt::T) where T
    t = range(t_span..., step = Δt)  # 時間軸
    x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

    x[1] = x₀  # 初期値
    for i in tqdm(1:length(x)-1)
        k₁ = dx(x[i])
        k₂ = dx(x[i]+k₁*Δt/2)
        k₃ = dx(x[i]+k₂*Δt/2)
        k₄ = dx(x[i]+k₃*Δt)
        x[i+1] = x[i] .+ (k₁ .+ 2k₂ .+ 2k₃ .+k₄) .* Δt/6
    end

    t, x
end


"""シミュレーションのデータを格納"""
@with_kw struct Data{T}
    t::Vector{T}
    q::Vector{Vector{T}}
    q_dot::Vector{Vector{T}}
    q_dot_dot::Vector{Vector{T}}
    qd::Vector{Vector{T}}
    qd_dot::Vector{Vector{T}}
    qd_dot_dot::Vector{Vector{T}}
    τ::Vector{Vector{T}}
    error::Vector{T}
end

const R = 0.01

"""適当な目標アクチュエータ変位"""
function calc_qd(t::T) where T
    [
        R * sin(3t),
        R * sin(3t + π/2),
        R * sin(3t + π),
    ]
end


"""適当な目標アクチュエータ速度"""
function calc_qd_dot(t::T) where T
    [
        R * 3cos(3t),
        R * 3cos(3t + π/2),
        R * 3cos(3t + π),
    ]
end


"""適当な目標アクチュエータ加速度"""
function calc_qd_dot_dot(t::T) where T
    [
        R * -9sin(3t),
        R * -9sin(3t + π/2),
        R * -9sin(3t + π),
    ]
end


"""状態方程式"""
function X_dot(X::Vector{T}) where T
    q = X[1:3]
    q_dot = X[4:6]

    τ = [10.0, 0.0, 0.0]

    x1_dot = q_dot
    x2_dot = inv(M(q)) * (τ .- G(q) .- (C(q, q_dot) .+ D)*q_dot .- K*q)
    #x2_dot = M(q) \ (τ .- G(q) .- (C(q, q_dot) .+ D)*q_dot .- K*q)

    #println(norm([x1_dot; x2_dot]))
    return [x1_dot; x2_dot]

end


"""シミュレーション実行"""
function run_simulation(TIME_SPAN::T=1.0, TIME_INTERVAL::T=0.0001) where T


    # 制御器のパラメータ
    kinematic = KinematicController(Dynamics.K)
    PDandFL = PDandFBController(
        Matrix{Float64}(I, 3, 3)*200, Matrix{Float64}(I, 3, 3)*1e4
    )
    passivity = PassiveController(
        Matrix{Float64}(I, 3, 3)*100, [10.0, 10.0, 10.0],
    )


    P = kinematic


    q₀ = zeros(T, 3)
    q_dot₀ = zeros(T, 3)

    t = Vector(0.0:TIME_INTERVAL:TIME_SPAN)

    data = Data(
        t = t,
        q = Vector{typeof(q₀)}(undef, length(t)),
        q_dot = Vector{typeof(q₀)}(undef, length(t)),
        q_dot_dot = Vector{typeof(q₀)}(undef, length(t)),
        qd = Vector{typeof(q₀)}(undef, length(t)),
        qd_dot = Vector{typeof(q₀)}(undef, length(t)),
        qd_dot_dot = Vector{typeof(q₀)}(undef, length(t)),
        τ = Vector{typeof(q₀)}(undef, length(t)),
        error = Vector{T}(undef, length(t)),
    )


end



@time run_simulation()