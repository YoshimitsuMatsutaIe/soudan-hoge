"""制御器のまとめ
"""
module Controller

using LinearAlgebra

include("dynamics.jl")
using .Dynamics: M, C, G, K, D


export KinematicController
export PDandFBController
export PassiveController
export calc_torque


"""運動学制御器"""
struct KinematicController{T}
    K_kin::Matrix{T}
end


function calc_torque(p::KinematicController{T}, q::Vector{T}, qd::Vector{T}) where T
    ek = qd .- q
    return p.K_kin * ek
end



"""PD + フィードバック線形化制御器"""
struct PDandFBController{T}
    Kd::Matrix{T}
    Kp::Matrix{T}
end


function calc_torque(
    p::PDandFBController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T}
    ) where T

    α = M(q)
    β = C(q, q_dot)*q_dot .+ D*q_dot .+ K*q .+ G(q)

    q̃ = q .- qd
    q̃_dot = q_dot .- qd_dot

    τ_dash = qd_dot_dot .- p.Kd*q̃_dot .- p.Kp*q̃

    return α*τ_dash .+ β
end


"""受動性に基づく制御器"""
struct PassiveController{T}
    Λ::Matrix{T}
    KG::Matrix{T}
end



function calc_torque(
    p::PassiveController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T}
    ) where T

    q̃ = q .- qd
    q̃_dot = q_dot .- qd_dot

    v = qd_dot .- p.Λ*q̃
    a = qd_dot_dot .- p.Λ*q̃_dot
    r = q̃_dot .+ p.Λ*q̃
    
    return M(q)*a .+ C(q, q_dot)*v .+ G(q) .+ K*q .+ D*v .- p.KG*r
end


"""受動性に基づく適応制御"""
struct PassiveAdaptiveController{T}
    Γ::Matrix{T}
    KG::Matrix{T}
end



function calc_torque(
    p::PassiveAdaptiveController{T},
    q::Vector{T}, q_dot::Vector{T},
    qd::Vector{T}, qd_dot::Vector{T}, qd_dot_dot::Vector{T}
    ) where T

    return 


end



end