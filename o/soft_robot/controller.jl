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
using .Dynamics: M, C, G, K, D, uncertain_K, uncertain_D, invM


export KinematicController
export PDandFBController
export PassivityBasedController
export PassivityBasedAdaptiveController
export SDREController
export MPCController
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
    n::T  # 予測ホライズン
    Δt::T  # 予測ホライズンの刻み時価
    isUncertainty::Bool
end


"""線形化したときのA行列"""
function A!(
    q::Vector{Float64}, q_dot::Vector{Float64}, H::Vector{Float64},
    out::Matrix{Float64}
    )
    ccall(
        (:q_dot_dot, "o/soft_robot/derived/ikko_dake/eqs/c_so/A.so"),
        Cvoid,
        (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
        H[1], H[2], H[3], q[1], q_dot[1], q[2], q_dot[2], q[3], q_dot[3], out
    )
end


"""線形化したときのA行列"""
function A(
    q::Vector{Float64}, q_dot::Vector{Float64}, H::Vector{Float64},
    )
    Z = Matrix{Float64}(undef, 9, 9)
    A!(q, q_dot, H, Z)
    Z
end


"""状態方程式のドリフト項"""
function fx!(
    q::Vector{Float64}, q_dot::Vector{Float64}, H::Vector{Float64},
    out::Vector{Float64}
    )
    ccall(
        (:q_dot_dot, "o/soft_robot/derived/ikko_dake/eqs/c_so/fx.so"),
        Cvoid,
        (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
        H[1], H[2], H[3], q[1], q_dot[1], q[2], q_dot[2], q[3], q_dot[3], out
    )
end


"""状態方程式のドリフト項"""
function fx(
    q::Vector{Float64}, q_dot::Vector{Float64}, H::Vector{Float64},
    )
    Z = Vector{Float64}(undef, 9)
    fx!(q, q_dot, H, Z)
    Z
end


function ℱ(A::Matrix{T}, n::Int64) where T
    m = size(A, 1)  # システムの次元
    Z = Marix{T}(undef, m*n, m)

    for i in 1:n
        Z[(i-1)*m+1:i*m, :] = A^i
    end
    
    Z
end


function 𝒢(A::Matrix{T}, B::Matrix{T}, n::Int64) where T
    m = size(A, 1)  # システムの次元
    _m = size(A, 2)
    Z = Marix{T}(undef, m*n, _m*n)

    for i in 1:n
        for j in 1:n
            if j > i
                Z[(i-1)*m+1:i*m, (j-1)*_m+1:j*_m] = zero(B)
            else
                Z[(i-1)*m+1:i*m, (j-1)*_m+1:j*_m] = A^(i-1) * B
            end
        end
    end
    
    Z
end


function 𝒮(A::Matrix{T}, n::Int64) where T
    m = size(A, 1)  # システムの次元
    Z = Marix{T}(undef, m*n, m*n)

    for i in 1:n
        for j in 1:n
            if j > i
                Z[(i-1)*m+1:i*m, (j-1)*m+1:j*m] = zero(A)
            else
                Z[(i-1)*m+1:i*m, (j-1)*m+1:j*m] = A^(i-1)
            end
        end
    end
    
    Z
end


function ℋ(C::Matrix{T}, n::Int64) where T
    m = size(C, 1)
    Z = zeros(T, m, m)
    for i in 1:n
        Z[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = C
    end

    Z
end


"""MPCで入力を計算"""
function calc_torque(
    p::MPCController{T},
    q::Vector{T}, q_dot::Vector{T}, H::Vector{T},
    qd, qd_dot,
    t::T
    ) where T

    X₀ = [q; q_dot; H]

    # A, B，C行列を計算
    A = A(q, q_dot, H)
    B = [
        zeros(T, 3, 3)
        invM(q)
        zeros(T, 3, 3)
    ]
    C = Matrix{T}(I, 9, 9)

    # 目標状態ベクトルと目標入力ベクトル作成
    Yref = Vector{T}(undef, 9*p.n)
    for i in 1:p.n
        Yref[(i-1)*9+1:i*9] = [
            qd(t + i*p.Δt)
            qd_dot(t + i*p.Δt)
            zeros(T, 3)
        ]
    end
    Uref = zero(Yref)

    # 外乱ベクトル作成
    W = Vector{T}(undef, 9*p.n)
    for i in 1:p.n
        W[(i-1)*9+1:i*9] = [
            zeros(T, 6)
            H
        ]
    end

    # F, G, S行列作成
    ℱ = ℱ(A, p.n)
    𝒢 = 𝒢(A, B, p.n)
    𝒮 = 𝒮(A, p.n)

    ℳ = 𝒢' * C' * p.Q * C * G .+ R
    𝒩 = (C*(ℱ*X₀ .+ 𝒮*W .- Yref))' * Q * C * 𝒢 .- Uref'*R
    Uopt = -inv(ℳ) * 𝒩'  # 最適入力

    return Uopt[1:3]
end




end