using Plots
using LinearAlgebra
using ProgressBars
using StaticArrays
using ArraysOfArrays


function split_vec_of_arrays(u)
    vec.(u) |>
    x ->
    VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end

function solve_RungeKutta(dx, x₀::Vector{T}, t_span, Δt::T) where T
    """ルンゲクッタ法（4次）"""

    t = range(t_span..., step = Δt)  # 時間軸
    x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

    x[1] = x₀  # 初期値
    for i in tqdm(1:length(x)-1)
        k₁ = dx(x[i])
        k₂ = dx(x[i]+k₁*Δt/2)
        k₃ = dx(x[i]+k₂*Δt/2)
        k₄ = dx(x[i]+k₃*Δt)
        x[i+1] = x[i] .+ (k₁ .+ 2k₂ .+ 2k₃ .+k₄) .* Δt/6
    end

    t, x
end


const N = 1  # セクションの数


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
    Z = Matrix{Float64}(undef, 3*N, 3*N)
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
    Z = Matrix{Float64}(undef, 3*N, 3*N)

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
    Z = Vector{Float64}(undef, 3*N)

    G!(q, Z)

    return Z
end

const K = diagm([1700.0, 1700.0, 1700.0])
const D = diagm([110.0, 110.0, 110.0])

"""状態方程式"""
function X_dot(X::Vector{T}) where T
    q = X[1:3]
    q_dot = X[4:6]
    xi = 1.0

    τ = [0.0, 0.0, 0.0]

    x1_dot = q_dot
    x2_dot = inv(M(q)) * (τ .- G(q) .- (C(q, q_dot) .+ D)*q_dot .- K*q)
    #x2_dot = M(q) \ (τ .- G(q) .- (C(q, q_dot) .+ D)*q_dot .- K*q)

    #println(norm([x1_dot; x2_dot]))
    return [x1_dot; x2_dot]

end


"""合ってるかテスト"""
function test()
    println("計算中...")
    TIME_SPAN = 5.0
    TIME_INTERVAL = 0.0001  # これより大きいと発散．方程式が硬すぎる?

    q = [0.0, 0.0, 0.0]
    xi = 1.0
    q_dot = [0.0, 0.0, 0.0]

    t, X = solve_RungeKutta(
        X_dot,
        [q; q_dot],
        (0.0, TIME_SPAN),
        TIME_INTERVAL
    )

    l1, l2, l3, l1_dot, l2_dot, l3_dot = split_vec_of_arrays(X)
    fig = plot(xlim=(0, TIME_SPAN))
    plot!(fig, t, l1, label="l1")
    plot!(fig, t, l2, label="l2")
    plot!(fig, t, l3, label="l3")
    

    fig2 = plot(xlim=(0, TIME_SPAN))
    plot!(fig2, t, l1_dot, label="l1_dot")
    plot!(fig2, t, l2_dot, label="l2_dot")
    plot!(fig2, t, l3_dot, label="l3_dot")

    fig_I = plot(fig, fig2, layout=(1, 2), size=(800, 400))
    savefig(fig_I, "julia.png")



    # アニメ制作



end


@time test()