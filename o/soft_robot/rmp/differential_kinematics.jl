"""微分運動学

"""
module DifferentialKinematics


export J_0
export J_1
export J_2
export J_3
export J_4

const N = 5


function J_0!(q::Vector{Float64}, ξ::Float64, out::Matrix{Float64})
    ccall(
        (:J_0, "o/soft_robot/derivation_of_kinematics/derived/so/J_0.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_0(q::Vector{Float64}, ξ::Float64, n::Int64)
    Z = Matrix{Float64}(undef, 3*N, 3)
    J_0!(q, ξ, Z)
    return Z'[:, 1:3*n]
end


function J_1!(q::Vector{Float64}, ξ::Float64, out::Matrix{Float64})
    ccall(
        (:J_1, "o/soft_robot/derivation_of_kinematics/derived/so/J_1.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_1(q::Vector{Float64}, ξ::Float64, n::Int64)
    Z = Matrix{Float64}(undef, 3*N, 3)
    J_1!(q, ξ, Z)
    return Z'[:, 1:3*n]
end


function J_2!(q::Vector{Float64}, ξ::Float64, out::Matrix{Float64})
    ccall(
        (:J_2, "o/soft_robot/derivation_of_kinematics/derived/so/J_2.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_2(q::Vector{Float64}, ξ::Float64, n::Int64)
    Z = Matrix{Float64}(undef, 3*N, 3)
    J_2!(q, ξ, Z)
    return Z'[:, 1:3*n]
end



function J_3!(q::Vector{Float64}, ξ::Float64, out::Matrix{Float64})
    ccall(
        (:J_3, "o/soft_robot/derivation_of_kinematics/derived/so/J_3.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_3(q::Vector{Float64}, ξ::Float64, n::Int64)
    Z = Matrix{Float64}(undef, 3*N, 3)
    J_3!(q, ξ, Z)
    return Z'[:, 1:3*n]
end


function J_4!(q::Vector{Float64}, ξ::Float64, out::Matrix{Float64})
    ccall(
        (:J_4, "o/soft_robot/derivation_of_kinematics/derived/so/J_4.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_4(q::Vector{Float64}, ξ::Float64, n::Int64)
    Z = Matrix{Float64}(undef, 3*N, 3)
    J_4!(q, ξ, Z)
    return Z'[:, 1:3*n]
end










end


using .DifferentialKinematics


using LinearAlgebra


q = [
    0, 0, 0,
    0, 0, 0,
    0.1, 0.2, 0.1
]

q0 = [0.01, 0.05, 0.1]

J = J_2(q, 1.0, 3)
J0 = J_0(q0, 1.0, 1)