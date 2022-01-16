"""制御器のまとめ
"""
module Controller

using LinearAlgebra
using Parameters  # 構造体にキーワード引数をつけるためのモジュール
using MatrixEquations  #　行列方程式ソルバー
using ControlSystems  # 制御関係
using Optim  # 最適化ライブラリ
using JuMP

include("dynamics.jl")
using .Dynamics: M, C, G, K, D, uncertain_K, uncertain_D


export KinematicController
export PDandFBController
export PassivityBasedController
export PassivityBasedAdaptiveController
export SDREController
export calc_torque
export θp_dot



"""運動学制御器"""
@with_kw struct KinematicController{T}
    K_kin::Matrix{T}
    isUncertainty::Bool
end


"""入力を計算"""
function calc_torque(
    p::KinematicController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T},
    ) where T
    ek = qd .- q
    return p.K_kin * ek
end





"""PD + フィードバック線形化制御器"""
@with_kw struct PDandFBController{T}
    Kd::Matrix{T}
    Kp::Matrix{T}
    isUncertainty::Bool
end


"""入力を計算"""
function calc_torque(
    p::PDandFBController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T},
    ) where T
    
    α = M(q)
    q̃ = q .- qd
    q̃_dot = q_dot .- qd_dot

    τ_dash = qd_dot_dot .- p.Kd*q̃_dot .- p.Kp*q̃
    
    if p.isUncertainty  # 不確かさあり
        β = C(q, q_dot)*q_dot .+ uncertain_D*q_dot .+ uncertain_K*q .+ G(q)
    else
        β = C(q, q_dot)*q_dot .+ D*q_dot .+ K*q .+ G(q)
    end

    α*τ_dash .+ β
end


"""受動性に基づく制御器"""
@with_kw struct PassivityBasedController{T}
    Λ::Matrix{T}
    KG::Matrix{T}
    isUncertainty::Bool
end


"""入力を計算"""
function calc_torque(
    p::PassivityBasedController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T},
    ) where T
    #print("hoge")
    q̃ = q .- qd
    q̃_dot = q_dot .- qd_dot

    v = qd_dot .- p.Λ*q̃
    a = qd_dot_dot .- p.Λ*q̃_dot
    r = q̃_dot .+ p.Λ*q̃
    
    if p.isUncertainty  # 不確かさあり
        return M(q)*a .+ C(q, q_dot)*v .+ G(q) .+ uncertain_K*q .+ uncertain_D*v .- p.KG*r
    else
        return M(q)*a .+ C(q, q_dot)*v .+ G(q) .+ K*q .+ D*v .- p.KG*r
    end
end


"""受動性に基づく適応制御"""
@with_kw struct PassivityBasedAdaptiveController{T}
    invΓ::Matrix{T}
    Λ::Matrix{T}
    KG::Matrix{T}
    isUncertainty::Bool
end


"""パラメータの更新式"""
function θp_dot(
    p::PassivityBasedAdaptiveController{T},
    q::Vector{T}, q_dot::Vector{T}, qd::Vector{T}, qd_dot::Vector{T},
    ) where T

    q̃ = q .- qd
    q̃_dot = q_dot .- qd_dot

    v = qd_dot .- p.Λ*q̃
    r = q̃_dot .+ p.Λ*q̃

    Y = [diagm(q) diagm(v)]

    return -p.invΓ * Y' * r
end


"""入力を計算"""
function calc_torque(
    p::PassivityBasedAdaptiveController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T},
    θp::Vector{T}
    ) where T
    #print("hoge")

    q̃ = q .- qd
    q̃_dot = q_dot .- qd_dot

    v = qd_dot .- p.Λ*q̃

    r = q̃_dot .+ p.Λ*q̃

    Y = [diagm(q) diagm(v)]
    

    return Y * θp .- p.KG*r


end



"""SDRE法"""
@with_kw struct SDREController{T}
    Q::Matrix{T}
    R::Matrix{T}
    isUncertainty::Bool
end


"""可制御性判定"""
function isControllable(A, B)
    n, _ = size(A)
    _, m = size(B)
    #println(n, m)
    Co = Matrix{T}(undef, n, n*m)
    for i in 1:n
        Co[:, m*(i-1)+1:m*i] = (A^(i-1))*B
    end
    return rank(Co) == n
end


"""SDRE法の入力を計算"""
function calc_torque(
    p::SDREController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T},
    ) where T
    #println(uncertain_K)
    x = [
        qd .- q
        qd_dot .- q_dot
    ]
    invM = inv(M(q))

    if p.isUncertainty  # 不確かさありのとき
        A = [
            zeros(T, 3, 3) Matrix{T}(I, 3, 3)
            -invM*uncertain_K -invM*(C(q, q_dot) .+ uncertain_D)
        ]
        B = [
            zeros(T, 3, 3)
            invM
        ]
    else  # 不確かさ無し
        #println(K)
        A = [
            zeros(T, 3, 3) Matrix{T}(I, 3, 3)
            -invM*K -invM*(C(q, q_dot) .+ D)
        ]
        B = [
            zeros(T, 3, 3)
            invM
        ]
    end

    #println("Co = \n", ctrb(A, B))
    #println("rank = ", ctrb(A, B) |> rank)

    if rank(ctrb(A, B)) == size(A, 1)
        println("可制御!")
    else
        println("ランク落ち...")
    end

    P, _, _ = arec(A, B, p.R, p.Q, zero(B))  # リカッチ方程式を解く
    #println(P)
    opt_gain = inv(p.R) * B' * P  # 最適フィードバックゲイン
    println(-opt_gain*x .- G(q))
    return -opt_gain*x .- G(q)
end


"""MPC制御器"""
@with_kw struct MPCController{T}
    Q::Matrix{T}
    R::Matrix{T}
    predit_span::T  # 予測ホライズン
    isUncertainty::Bool
end


function state_eq!(X_dot::Vector{T}, X::Vector{T}, p, t::T) where T
    q = X[1:3]
    q_dot = X[4:6]
    H = X[7:9]
    #println(t)
    τ = p
    X_dot[1:3] = q_dot
    X_dot[4:6] = Dynamics.calc_q_dot_dot(τ, q, q_dot, H)
    X_dot[7:9] = H_dot(H, q, q_dot)
end



"""MPCのコスト関数"""
function cost(
    p::MPCController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T},
    ) where T



end

function calc_torque(
    p::MPCController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T},
    ) where T




end




end