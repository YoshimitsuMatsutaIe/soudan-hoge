"""運動学

"""
module Kinematics

export Phi_0
export Phi_1
export Phi_2

export Arm
ctrlab2021_soudan\o\soft_robot\derivation_of_kinematics\derived\N_is_3\c_src\J_s\J_0.so

function Phi_0!(q::Vector{Float64}, ξ::Vector{Float64}, out::Vector{Float64})
    ccall(
        (:Phi_0, "/o/soft_robot/derivation_of_kinematics/derived/N_is_3/c_src/Phi_s/Phi_0.so"),
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
        (:Phi_1, "/o/soft_robot/derivation_of_kinematics/derived/N_is_3/c_src/Phi_s/Phi_1.so"),
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
        (:Phi_2, "/o/soft_robot/derivation_of_kinematics/derived/N_is_3/c_src/Phi_s/Phi_2.so"),
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

"""アームのディスク位置を計算"""
function Arm(q::Vector{T}, m::Int64) where T
    n = length(Ξ)
    Phis = Vector{Vector{T}}(undef, n)
    for (i, ξ) in enumerate(Ξ)
        if m == 0
            Phis[i] = Phi_0(q, ξ)
        elseif m == 1
            Phis[i] = Phi_1(q, ξ)
        else
            Phis[i] = Phi_2(q, ξ)
        end


    end
    return Phis
end

end

# using .Kinematics

# q = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# p0 = Phi_0(q, 0.5)
# p1 = Phi_1(q, 0.5)
# p2 = Phi_2(q, 0.5)

# Arm(q, 0)
# println(p0)
# println(p1)
# println(p2)