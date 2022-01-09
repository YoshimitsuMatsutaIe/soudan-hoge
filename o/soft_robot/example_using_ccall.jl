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

"""ルンゲクッタ法（4次）"""
function solve_RungeKutta(dx, x₀::Vector{T}, t_span, Δt::T) where T
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

    τ = [10.0, 0.0, 0.0]

    x1_dot = q_dot
    x2_dot = inv(M(q)) * (τ .- G(q) .- (C(q, q_dot) .+ D)*q_dot .- K*q)
    #x2_dot = M(q) \ (τ .- G(q) .- (C(q, q_dot) .+ D)*q_dot .- K*q)

    #println(norm([x1_dot; x2_dot]))
    return [x1_dot; x2_dot]

end


const Ξ = Vector(0:0.01:1)


"""1フレームを描写"""
function draw_frame(t::T, q::Vector{T}, q_dot::Vector{T}, fig_shape) where T
    

    arm = Arm(q, Ξ)

    x, y, z = split_vec_of_arrays(arm)
    fig = plot(
        x, y, z,
        #marker=:circle,
        aspect_ratio = 1,
        #markersize=2,
        label="arm",
        xlabel = "X[m]", ylabel = "Y[m]", zlabel = "Z[m]",
    )

    # scatter!(
    #     fig,
    #     [xd_true[1]], [xd_true[2]],
    #     label="xd_true",
    #     markershape=:star6,
    # )

    # scatter!(
    #     fig,
    #     [xd[1]], [xd[2]],
    #     label="xd_true",
    #     markershape=:star6,
    # )


    plot!(
        fig,
        xlims=(fig_shape.xl, fig_shape.xu),
        ylims=(fig_shape.yl, fig_shape.yu),
        zlims=(fig_shape.zl, fig_shape.zu),
        legend = true,
        size=(600, 600),
        title = string(round(t, digits=2)) * "[s]"
    )

    return fig
end


"""アニメ作成"""
function make_animation(t, x)
    println("アニメ作成中...")
    # 枚数決める
    #println(data.t)
    epoch_max = 100
    epoch = length(t)
    if epoch < epoch_max
        step = 1
    else
        step = div(epoch, epoch_max)
    end

    #println(step)

    x_max = 0.18
    x_min = -0.18
    y_max = 0.18
    y_min = -0.18
    z_max = 0.15
    z_min = 0.0
    max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)*0.5
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2
    fig_shape = (
        xl = x_mid-max_range, xu = x_mid+max_range,
        yl = y_mid-max_range, yu = y_mid+max_range,
        zl = z_mid-max_range, zu = z_mid+max_range,
    )

    anim = Animation()
    @gif for i in tqdm(1:step:length(t))
        _fig = draw_frame(t[i], x[i][1:3], x[i][4:6], fig_shape)
        frame(anim, _fig)
    end



    gif(anim, "julia.gif", fps = 60)

    println("アニメ作成完了")
end


"""合ってるかテスト"""
function test()
    println("計算中...")
    TIME_SPAN = 10.0
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
    make_animation(t, X)


end

#@time X = Arm([0.0, 0.0, 0.0], Ξ)
@time test()