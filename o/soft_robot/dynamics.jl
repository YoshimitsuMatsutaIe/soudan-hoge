"""動力学関連
"""
module Dynamics

using LinearAlgebra

"""Mのラッパー"""
function M!(q::Vector{Float64}, out::Matrix{Float64})
    ccall(
        (:M, "o/soft_robot/derived/ikko_dake/eqs/c_so/M.so"),
        Cvoid,
        (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
        q[1], q[2], q[3], out
    )
end


"""慣性行列"""
function M(q::Vector{Float64})
    Z = Matrix{Float64}(undef, 3, 3)
    M!(q, Z)
    return Z
end


function C!(q::Vector{Float64}, q_dot::Vector{Float64}, out::Matrix{Float64})
    ccall(
        (:C, "o/soft_robot/derived/ikko_dake/eqs/c_so/C.so"),
        Cvoid,
        (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
        q[1], q[2], q[3], q_dot[1], q_dot[2], q_dot[3], out
    )
end


"""コリオリ＋遠心力行列"""
function C(q::Vector{Float64}, q_dot::Vector{Float64})
    Z = Matrix{Float64}(undef, 3, 3)

    C!(q, q_dot, Z)

    return Z
end


"""Gのラッパー"""
function G!(q::Vector{Float64}, out::Vector{Float64})
    ccall(
        (:G, "o/soft_robot/derived/ikko_dake/eqs/c_so/G.so"),
        Cvoid,
        (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
        q[1], q[2], q[3], out
    )
end

"""重力行列"""
function G(q::Vector{Float64})
    Z = Vector{Float64}(undef, 3)

    G!(q, Z)

    return Z
end


const K = diagm([1700.0, 1700.0, 1700.0])  # 剛性行列の真値
const D = diagm([110.0, 110.0, 110.0])  # 減衰行列の真値

const uncertain_K = diagm([1020.0, 1020.0, 1020.0])  # 不確かな剛性行列
const uncertain_D = diagm([77.0, 77.0, 77.0])  # 不確かな減衰行列

const αh = 23.705
const βh = 1.7267
const γh = -42.593

function h_dot(h::T, q::T, q_dot::T) where T
    q * (αh - (βh * sign(q_dot*h) + γh)*abs(h))
end


"""ヒステリシスベクトルの更新式"""
function H_dot(H::Vector{T}, q::Vector{T}, q_dot::Vector{T}) where T
    [
        h_dot(H[1], q[1], q_dot[1]),
        h_dot(H[2], q[2], q_dot[2]),
        h_dot(H[3], q[3], q_dot[3]),
    ]
end


"""加速度"""
function calc_q_dot_dot(
    τ::Vector{T}, q::Vector{T}, q_dot::Vector{T},
    H::Vector{T}, 
    ) where T
    inv(M(q)) * (τ .- (C(q, q_dot) .+ D)*q_dot .- K*q .- G(q) .- H) |> vec
end

end