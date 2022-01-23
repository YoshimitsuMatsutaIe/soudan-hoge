"""運動学

"""
module Kinematics

export Phi_0
export Phi_1
export Phi_2
export Phi_3
export Phi_4

export Arm
export Ξ

function Phi_0!(q::Vector{Float64}, ξ::Float64, out::Vector{Float64})
    ccall(
        (:Phi_0, "o/soft_robot/derivation_of_kinematics/derived/so/Phi_0.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_0(q::Vector{Float64}, ξ::Float64)
    Z = Vector{Float64}(undef, 3)
    Phi_0!(q, ξ, Z)
    return Z
end

function Phi_1!(q::Vector{Float64}, ξ::Float64, out::Vector{Float64})
    ccall(
        (:Phi_1, "o/soft_robot/derivation_of_kinematics/derived/so/Phi_1.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_1(q::Vector{Float64}, ξ::Float64)
    Z = Vector{Float64}(undef, 3)
    Phi_1!(q, ξ, Z)
    return Z
end

function Phi_2!(q::Vector{Float64}, ξ::Float64, out::Vector{Float64})
    ccall(
        (:Phi_2, "o/soft_robot/derivation_of_kinematics/derived/so/Phi_2.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_2(q::Vector{Float64}, ξ::Float64)
    Z = Vector{Float64}(undef, 3)
    Phi_2!(q, ξ, Z)
    return Z
end


function Phi_3!(q::Vector{Float64}, ξ::Float64, out::Vector{Float64})
    ccall(
        (:Phi_3, "o/soft_robot/derivation_of_kinematics/derived/so/Phi_3.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_3(q::Vector{Float64}, ξ::Float64)
    Z = Vector{Float64}(undef, 3)
    Phi_3!(q, ξ, Z)
    return Z
end

function Phi_4!(q::Vector{Float64}, ξ::Float64, out::Vector{Float64})
    ccall(
        (:Phi_4, "o/soft_robot/derivation_of_kinematics/derived/so/Phi_4.so"),
        Cvoid,
        (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
        q, ξ, out
    )
end

"""グローバル位置"""
function Phi_4(q::Vector{Float64}, ξ::Float64)
    Z = Vector{Float64}(undef, 3)
    Phi_4!(q, ξ, Z)
    return Z
end



const Ξ = Vector(0:0.01:1)


"""アームのディスク位置を計算"""
function Arm(q::Vector{T}, section_n::Int64) where T
    n = length(Ξ)
    Phis = Vector{Vector{T}}(undef, section_n*n)
    for (i, ξ) in enumerate(Ξ)
        Phis[i] = Phi_0(q, ξ)
        if section_n > 1
            Phis[i+n] = Phi_1(q, ξ)
        elseif section_n > 2
            Phis[i+2n] = Phi_2(q, ξ)
        elseif  section_n > 3
            Phis[i+3n] = Phi_3(q, ξ)
        elseif  section_n > 4
            Phis[i+4n] = Phi_4(q, ξ)
        end
        #println(Phis[i], Phis[i+n], Phis[i+2n])
    end
    return Phis
end

# """アームのディスク位置を計算"""
# function Arm(q::Vector{T}, m::Int64) where T
#     n = length(Ξ)
#     Phis = Vector{Vector{T}}(undef, n)
#     for (i, ξ) in enumerate(Ξ)
#         if m == 0
#             Phis[i] = Phi_0(q, ξ)
#         elseif m == 1
#             Phis[i] = Phi_1(q, ξ)
#         else
#             Phis[i] = Phi_2(q, ξ)
#         end


#     end
#     return Phis
# end

end

# using .Kinematics

# q = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# p0 = Phi_0(q, 0.5)
# p1 = Phi_1(q, 0.5)
# p2 = Phi_2(q, 0.5)

# #Arm(q, 0)
# println(p0)
# println(p1)
# println(p2)