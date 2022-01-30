


"""
運動学
"""
module Kinematics


export calc_local
export calc_global


const c1 = 837019575
const c2 = 4133430
const c3 = 32805
const c4 = 486
const c5 = 18
const c6 = 55801305
const c7 = 688905
const c8 = 3645
const c9 = 81

const c10 = 279006525
const c11 = 1377810
const c12 = 10935
const c13 = 162

const c14 = 243
const c15 = 2066715

const r = 0.0125
const L0 = 0.15

const sq3 = sqrt(3)



"""Aを計算"""
function As(q::T) where T

    A1 = q[1]^2 + q[2]^2 + q[3]^2 - q[1]*q[2] - q[1]*q[3] - q[2]*q[3]
    A2 = 2*q[1] - q[2] - q[3]
    A3 = q[2] - q[3]
    A4 = 3*L0 + q[1] + q[2] + q[3]
    
    A1, A2, A3, A4
end



function P1(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    -(A2 * A1^4 * A4 * ξ^10) / ((c1 * r^9)) +
    (A2 * A1^3 * A4 * ξ^8) / (c2 * r^7) -
    (A2 * A1^2 * A4 * ξ^6) / (c3 * r^5) +
    (A2 * A1 * A4 * ξ^4) / (c4 * r^3)
    (A2 * A4 * ξ^2) / (c5 * r)
end


function P2(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    -(sq3 * A4 * A3 * A1^4 * ξ^10) / (c1 * r^9) +
    (sq3 * A4 * A3 * A1^3 * ξ^8) / (c2 * r^7) -
    (sq3 * A4 * A3 * A1^2 * ξ^6) / (c3 * r^5) +
    (sq3 * A4 * A1 * A2 * ξ^4) / (c4 * r^3) -
    (sq3 * A4 * A3 * ξ^2) / (c5 * r)
end


function P3(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    (2 * A1^4 * A4 * ξ^9) / (c6 * r^8) -
    (4 * A1^3 * A4 * ξ^7) / (c7 * r^6) +
    (2 * A1^2 * A4 * ξ^5) / (c8 * r^4) -
    (2 * A1 *A4 * ξ^3) / (c9 * r^2) +
    (A4 * ξ) / 3
end



function P(q::T, ξ::U) where {T, U}
    A1, A2, A3, A4 = As(q)
    
    [
        P1(A1, A2, A3, A4, ξ)
        P2(A1, A2, A3, A4, ξ)
        P3(A1, A2, A3, A4, ξ)
    ]
end



function R11(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    1 - (A2^2 * A1^4 * ξ^10) / (c1 * r^10) +
    (A2^2 * A1^3 * ξ^8) / (c1 * r^8) -
    (A2^2 * A1^2 * ξ^6) / (c3 * r^6) +
    (A1 * A2^2 * ξ^4) / (c4 * r^4) -
    (A2^2 * ξ^2) / (c5 * r^2)
end

function R12(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    (sq3 * A2 * A3 * A1^4 * ξ^10) / (c1 * r^10) +
    (sq3 * A2 * A3 * A1^3 * ξ^8) / (c2 * r^8) -
    (sq3 * A2 * A3 * A1^2 * ξ^6) / (c3 * r^6) +
    (sq3 * A2 * A3 * A1 * ξ^4) / (c4 * r^4) -
    (sq3 * A2 * A3 * ξ^2) / (c5 * r^2)
end

function R13(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    -(2 * A2 * A1^4 * ξ^9) / (c6 * r^9) +
    (4 * A2 * A1^3 * ξ^7) / (c7 * r^7) -
    (2 * A2 * A1^2 * ξ^5) / (c8 * r^5) +
    (2 * A2 * A1 * ξ^3) / (c9 * r^3) -
    (A2 * ξ) / (3 * r)
end

function R21(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    R12(A1, A2, A3, A4, ξ)
end

function R22(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    1 - (A3^2 * A1^4 * ξ^10) / (c10 * r^10) +
    (A3^2 * A1^3 * ξ^8) / (c11 * r^8) -
    (A3^2 * A1^2 * ξ^6) / (c12 * r^6) +
    (A3^2 * A1 * ξ^4) / (c13 * r^4) -
    (A3^2 * ξ^2) / (6 * r^2)
end

function R23(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    -(2*sq3 * A3 * A1^4 * ξ^9) / (c6 * r^9) +
    (4*sq3 * A3 * A1^3 * ξ^7) / (c7 * r^7) -
    (2*sq3 * A3 * A1^2 * ξ^5) / (c8 * r^5) +
    (2*sq3 * A3 * A1 * ξ^3) / (c9 * r^3) -
    (sq3 * A3 * ξ) / (3 * r)
end

function R31(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    -R13(A1, A2, A3, A4, ξ)
end

function R32(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    -R23(A1, A2, A3, A4, ξ)
end

function R33(A1::T, A2::T, A3::T, A4::T, ξ::T) where T
    1 - (2 * ξ^2 * A1) / (9 * r^2) +
    (2 * ξ^4 * A1^2) / (c14 * r^4) -
    (4 * ξ^6 * A1^3) / (c3 * r^6) +
    (2 * ξ^8 * A1^4) / (c15 * r^8) -
    (4 * ξ^10 * A1^5) / (c1 * r^10)
end

function R(q::T, ξ::U) where {T, U}
    A1, A2, A3, A4 = As(q)

    [
        R11(A1, A2, A3, A4, ξ) R12(A1, A2, A3, A4, ξ) R13(A1, A2, A3, A4, ξ)
        R21(A1, A2, A3, A4, ξ) R22(A1, A2, A3, A4, ξ) R23(A1, A2, A3, A4, ξ)
        R31(A1, A2, A3, A4, ξ) R32(A1, A2, A3, A4, ξ) R33(A1, A2, A3, A4, ξ)
    ]
end


function calc_local(N::Int64, Q::T, Ξ::T) where T

    Ps = Vector{Vector{typeof(T)}}(undef, N)
    Rs = Vector{Matrix{typeof(T)}}(undef, N)
    for i in 1:N
        Ps[i] = P(Q[3(i-1)+1:3i], Ξ[i])
        println(R(Q[3(i-1)+1:3i], Ξ[i]))
        Rs[i] = R(Q[3(i-1)+1:3i], Ξ[i])
    end

    Ps, Rs
end



function calc_global(N::Int64, Ps::Vector{T}, Rs::Vector{U}) where {T, U}

    Φs = Vector{T}(undef, N)
    Θs = Vector{U}(undef, N)

    Φs[1] = Ps[1]
    Θs[1] = Rs[1]

    for i in 2:N
        Θs[i] = Θs[i-1] * Rs[i]
        Φs[i] = Φs[i-1] + Θs[i-1]*Ps[i]
    end

    Φs, Θs
end




end