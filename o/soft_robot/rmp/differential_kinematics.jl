"""微分運動学

"""
module DifferentialKinematics


export J_0
export J_1
export J_2
export J_3
#export J_4




function J_0!(q::Vector{Float64}, ξ::Vector{Float64}, out::Matrix{Float64})
    ccall(
        (:J_0, "o/soft_robot/derivation_of_kinematics/derived/so/J_0.so"),
        Cvoid,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_0(q::Vector{Float64}, ξ::Float64, N::Int64)
    Z = Matrix{Float64}(undef, 3, N)
    J_0!(q, [ξ, 0, 0], Z)
    return Z
end


function J_1!(q::Vector{Float64}, ξ::Vector{Float64}, out::Matrix{Float64})
    ccall(
        (:J_1, "o/soft_robot/derivation_of_kinematics/derived/so/J_1.so"),
        Cvoid,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_1(q::Vector{Float64}, ξ)
    Z = Matrix{Float64}(undef, 3, 9)
    J_1!(q, [0, ξ, 0], Z)
    return Z
end


function J_2!(q::Vector{Float64}, ξ::Vector{Float64}, out::Matrix{Float64})
    ccall(
        (:J_2, "o/soft_robot/derivation_of_kinematics/derived/so/J_2.so"),
        Cvoid,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_2(q::Vector{Float64}, ξ)
    Z = Matrix{Float64}(undef, 3, 9)
    J_2!(q, [0, 0, ξ], Z)
    return Z
end



function J_3!(q::Vector{Float64}, ξ::Vector{Float64}, out::Matrix{Float64})
    ccall(
        (:J_2, "o/soft_robot/derivation_of_kinematics/derived/so/J_3.so"),
        Cvoid,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""ヤコビ行列"""
function J_3(q::Vector{Float64}, ξ)
    Z = Matrix{Float64}(undef, 3, 9)
    J_2!(q, [0, 0, ξ], Z)
    return Z
end














end


using .DifferentialKinematics


using LinearAlgebra


q = [
    0, 0, 0,
    0, 0, 0,
    0.1, 0.2, 0.1
]


J = J_2(q, 1.0)