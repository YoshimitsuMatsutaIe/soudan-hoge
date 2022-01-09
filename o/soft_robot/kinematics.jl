"""運動学

"""
module Kinematics


"""Phi0のラッパー"""
function Phi0!(q::Vector{Float64}, ξ::Float64, out::Vector{Float64})
    ccall(
        (:Phi0, "o/soft_robot/derived/ikko_dake/eqs/c_so/Phi0.so"),
        Cvoid,
        (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
        q[1], q[2], q[3], ξ, out
    )
end

"""グローバル位置"""
function Phi0(q::Vector{Float64}, ξ)
    Z = Vector{Float64}(undef, 3*N)
    Phi0!(q, ξ, Z)
    return Z
end


"""アームのディスク位置を計算"""
function Arm(q::Vector{T}, Ξ::Vector{T}) where T
    Phis = Vector{typeof(q)}(undef, length(Ξ))
    for (i, ξ) in enumerate(Ξ)
        Phis[i] = Phi0(q, ξ)
    end
    return Phis
end



end