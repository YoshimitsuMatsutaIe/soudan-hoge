"""制御器のまとめ
"""
module Controller

using LinearAlgebra
using Parameters  # 構造体にキーワード引数をつけるためのモジュール
using MatrixEquations  #　行列方程式ソルバー


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
    
    if p.isUncertainty  # 不確かさあり
        β = C(q, q_dot)*q_dot .+ uncertain_D*q_dot .+ uncertain_K*q .+ G(q)
        τ_dash = qd_dot_dot .- p.Kd*q̃_dot .- p.Kp*q̃
        return α*τ_dash .+ β
    else
        β = C(q, q_dot)*q_dot .+ D*q_dot .+ K*q .+ G(q)
        τ_dash = qd_dot_dot .- p.Kd*q̃_dot .- p.Kp*q̃
        return α*τ_dash .+ β
    end
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


"""SDRE法の入力を計算"""
function calc_torque(
    p::SDREController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T},
    H::VectorT
    ) where T

    x = [
        qd .- q
        qs_dot .- q_dot
    ]
    invM = inv(M(q))

    if p.isUncertainty  # 不確かさありのとき
        println("hoge")
    else  # 不確かさ無し
        A = [
            zeros(T, 3, 3) Matrix{T}(I, 3, 3)
            -invM*K -insM*(C(q, q_dot) .+ D)
        ]
        B = [
            zeros(T, 3, 3)
            invM
        ]


    S = zero(B)  # 何か不明
    P, _, _ = arec(A, B, p.R, p.Q, S)  # リカッチ方程式を解く
    K = inv(p.R) * B' * P  # 最適フィードバックゲイン
    
    return -K*x .- G(q) .- H
end





end