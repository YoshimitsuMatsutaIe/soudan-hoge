"""運動学

"""
module Kinematics

export Phi_0
export Phi_1
export Phi_2

export Arm


function Phi_0!(q::Vector{Float64}, ξ::Vector{Float64}, out::Vector{Float64})
    ccall(
        (:Phi_0, "o/soft_robot/rmp/Phi_0.so"),
        Cvoid,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_0(q::Vector{Float64}, ξ)
    Z = Vector{Float64}(undef, 3)
    Phi_0!(q, [ξ, 0, 0], Z)
    return Z
end

function Phi_1!(q::Vector{Float64}, ξ::Vector{Float64}, out::Vector{Float64})
    ccall(
        (:Phi_1, "o/soft_robot/rmp/Phi_1.so"),
        Cvoid,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_1(q::Vector{Float64}, ξ)
    Z = Vector{Float64}(undef, 3)
    Phi_1!(q, [0, ξ, 0], Z)
    return Z
end

function Phi_2!(q::Vector{Float64}, ξ::Vector{Float64}, out::Vector{Float64})
    ccall(
        (:Phi_2, "o/soft_robot/rmp/Phi_2.so"),
        Cvoid,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_2(q::Vector{Float64}, ξ)
    Z = Vector{Float64}(undef, 3)
    Phi_2!(q, [0, 0, ξ], Z)
    return Z
end





const Ξ = Vector(0:0.01:1)


"""アームのディスク位置を計算"""
function Arm(q::Vector{T},) where T
    n = length(Ξ)
    Phis = Vector{Vector{T}}(undef, 3*n)
    for (i, ξ) in enumerate(Ξ)
        Phis[i] = Phi_0(q, ξ)
        Phis[i+n] = Phi_1(q, ξ)
        Phis[i+2n] = Phi_2(q, ξ)
        #println(Phis[i], Phis[i+n], Phis[i+2n])
    end
    return Phis
end



end

# using .Kinematics

# q = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# p0 = Phi_0(q, 0.5)
# p1 = Phi_1(q, 0.5)
# p2 = Phi_2(q, 0.5)

# println(p0)
# println(p1)
# println(p2)