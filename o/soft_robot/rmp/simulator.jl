"""シミュレーション"""


using LinearAlgebra  # 逆行列とノルムの計算に使用
using Plots  # グラフ作成
#using ProgressBars
#using Parameters
using DifferentialEquations  # 微分方程式のソルバ
using LaTeXStrings
using ProgressBars

include("utils.jl")
include("kinematics.jl")
include("differential_kinematics.jl")
#include("controller.jl")
#include("plot_utils.jl")

using .Kinematics
using .DifferentialKinematics
#using .PlotTool
#using .Controller

const sec_n = 3

const R = 0.15
const ω = 0.1
const offset_z = 0.45

"""適当な目標位置変位"""
function calc_xd(t::T) where T
    # [
    #     R * sin(ω*t),
    #     R * cos(ω*t),
    #     offset_z,
    # ]
    [
        0.2cos(ω*2*t)*cos(ω*t)
        0.2cos(ω*2*t)*sin(ω*t)
        offset_z
    ]
    #Phi_2([0.01, 0.01, 0, 0, 0, 0, 0.0, 0.02, 0.01], 1.0)
end


"""適当な目標位置速度"""
function calc_xd_dot(t::T) where T
    [
        R * 3cos(3t),
        R * -3sin(3t),
        0,
    ]
end


"""適当な目標位置加速度"""
function calc_xd_dot_dot(t::T) where T
    [
        R * -9sin(3t),
        R * -9cos(3t),
        0,
    ]
end

"""ソフトマックス関数"""
function soft_max(s::T, α::T) where T
    s + 1/α * log(1 + exp(-2 * α * s))
end

"""ソフト正規化関数"""
function soft_normal(v, alpha)
    return v ./ soft_max(norm(v), alpha)
end

"""空間を一方向に伸ばす計量"""
function metric_stretch(v, alpha)
    xi = soft_normal(v, alpha)
    return xi * xi'
end

"""基本の計量"""
function basic_metric_H(f::Vector{T}, alpha::T, beta::T) where T
    dim = length(f)
    return beta .* metric_stretch(f, alpha) + (1 - beta) .* Matrix{T}(I, dim, dim)
end




function pullbacked_rmp(f, M, J)
    pulled_f = J' * f
    pulled_M = J' * M * J
    return pulled_f, pulled_M
end

function state_eq!(
    X_dot::Vector{T}, X::Vector{T},
    p,
    t::T
    ) where T
    X_dot[1:9] = X[10:18]
    #println("t = ", t)
    J = J_2(X[1:9], 1.0, sec_n)
    
    #println(pinv(J))
    #println(J*X[10:18])
    # A = pinv(J) * calc_xd_dot_dot(t)
    # B = pinv(J) * (-p.D*(J*X[10:18] .- calc_xd_dot(t)))
    # C = pinv(J) * (-p.K*(Phi_2(X[1:9], 1.0) - calc_xd(t)))
    # X_dot[10:18] =  A .+ B .+ C
    # println("A = ", norm(A))
    # println("B = ", norm(B))
    #println("C = ", norm(C))

    # 目標制御
    gain = 900.0
    max_speed = 300.0
    sigma_H = 1.0
    sigma_W = 1.0
    damp_r = 0.01
    ddq_damp_r = 0.05

    z0 = calc_xd(t)
    z = Phi_2(X[1:9], 1.0)
    dz = J * X[10:18]
    damp = gain / max_speed
    a = gain .* soft_normal(z0.-z, damp_r) .- damp*dz

    dis = norm(z0 .- z)
    weight = exp(-dis ./ sigma_W)
    beta = 1.0 - exp(-1/2 * (dis / sigma_H)^2)
    M = weight .* basic_metric_H(a, ddq_damp_r, beta)

    pulled_f, pulled_M = pullbacked_rmp(a, M, J)
    root_f = pulled_f
    root_M = pulled_M


    # ジョイント制限回避
    gamma_p = 0.05
    gamma_d = 1.0
    lambda = 0.5
    q = X[1:9]
    dq = X[10:18]

    q_max = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    q_min = -[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    diags = (q_max .- q_min) .* (exp.(-q) ./ (1.0 .+ exp.(-q)).^2)
    D_sigma =  diagm(diags)

    z = gamma_p .* (-q) .- gamma_d .* dq
    a = inv(D_sigma) * z
    M = lambda .* Matrix{T}(I, 9, 9)
    
    @. root_f += a
    @. root_M += M
    
    ddq = pinv(pulled_M) * pulled_f

    X_dot[10:18] = ddq


end






"""1フレームを描写"""
function draw_frame(
    t::T, q::Vector{T}, xd::Union{Vector{T}, Nothing},
    fig_shape::Union{
        NamedTuple{(:xl, :xu, :yl, :yu, :zl, :zu), Tuple{T, T, T, T, T, T}},
        Nothing
    }
    ) where T
    

    fig = plot()
    arm = Vector{Vector{T}}(undef, length(Ξ))
    for (i, ξ) in enumerate(Ξ)
        arm[i] = Phi_0(q, ξ)
    end
    x, y, z = split_vec_of_arrays(arm)
    plot!(
        fig,
        x, y, z,
        #marker=:circle,
        aspect_ratio = 1,
        #markersize=2,
        label="1",
        xlabel = "X[m]", ylabel = "Y[m]", zlabel = "Z[m]",
    )

    arm = Vector{Vector{T}}(undef, length(Ξ))
    for (i, ξ) in enumerate(Ξ)
        arm[i] = Phi_1(q, ξ)
    end
    x, y, z = split_vec_of_arrays(arm)
    plot!(
        fig,
        x, y, z,
        #marker=:circle,
        aspect_ratio = 1,
        #markersize=2,
        label="2",
        xlabel = "X[m]", ylabel = "Y[m]", zlabel = "Z[m]",
    )

    arm = Vector{Vector{T}}(undef, length(Ξ))
    for (i, ξ) in enumerate(Ξ)
        arm[i] = Phi_2(q, ξ)
    end
    x, y, z = split_vec_of_arrays(arm)
    plot!(
        fig,
        x, y, z,
        #marker=:circle,
        aspect_ratio = 1,
        #markersize=2,
        label="3",
        xlabel = "X[m]", ylabel = "Y[m]", zlabel = "Z[m]",
    )



    if !isnothing(xd)
        scatter!(
            fig,
            [xd[1]], [xd[2]], [xd[3]],
            label="xd",
            markershape=:star6,
        )
    end

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


# fig = draw_frame(
#     0.0,
#     [
#         0.0, 0.0, 0.0,
#         0.0, 0.0, 0.00,
#         0.0, 0.15, 0.1
#     ], nothing, nothing)

"""アニメ作成"""
function make_animation(sol, xd_func)
    println("アニメ作成中...")
    # 枚数決める
    #println(data.t)
    epoch_max = 100
    epoch = length(sol.t)
    if epoch < epoch_max
        step = 1
    else
        step = div(epoch, epoch_max)
    end

    #println(step)

    x_max = 0.2
    x_min = -0.2
    y_max = 0.2
    y_min = -0.2
    z_max = 0.55
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

    xd_all = [xd_func(t) for t in sol.t]
    x, y, z = split_vec_of_arrays(xd_all)

    anim = Animation()
    @gif for i in tqdm(1:step:length(sol.t))
        _fig = draw_frame(sol.t[i], sol.u[i], xd_func(sol.t[i]), fig_shape)
        plot!(_fig, x, y, z, label="xd line",)
        frame(anim, _fig)
    end



    gif(anim, "julia.gif", fps = 60)

    println("アニメ作成完了")
end

"""微分方程式を解く"""
function run_simulation(
    TIME_SPAN::T
    ) where T

    #X₀ = zeros(T, 18)
    X₀ = [
        0, 0, 0,
        0, 0, 0,
        0.01, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
    ]
    t_span = (0.0, TIME_SPAN)
    p = (
        K = Matrix{T}(I, 3, 3)*0.1,
        D = Matrix{T}(I, 3, 3)*1,
    )
    prob = ODEProblem(state_eq!, X₀, t_span, p)
    solve(prob)#, tstops=[10000])
end


"""実行"""
function exmample()

    TIME_SPAN = 80.0

    fig0 = plot(xlims=(0.0, TIME_SPAN),) #ylims=(0,0.02))

    fig1 = plot(xlims=(0.0, TIME_SPAN))
    fig2 = plot(xlims=(0.0, TIME_SPAN))

    


    @time sol = run_simulation(TIME_SPAN)

    error = Vector{Vector{Float64}}(undef, length(sol.t))
    for (i, t) in enumerate(sol.t)
        error[i] = calc_xd(t) .- Phi_2(sol.u[i][1:9], 1.)
    end
    L2_error = norm.(error)
    plot!(
        fig0,
        sol.t, L2_error,
        label="err",
        legend=:outerright,
        #c=param.color,
        ylabel="error [m]",
        ylims=(0.0, maximum(L2_error))
    )

    plot!(
        fig1,
        sol, vars=[(0,1), (0,2), (0,3), (0,4), (0,5), (0,6), (0,7), (0,8), (0,9)],
        label=["l1_1" "l1_2" "l1_3" "l2_1" "l2_2" "l2_3" "l3_1" "l3_2" "l3_3"],
        legend=:outerright,
        #c=param.color,
        #linestyle=[:solid :dash :dot],
        ylabel="l [m]"
    )
    plot!(
        fig2,
        sol, vars=[(0,10), (0,11), (0,12), (0,13), (0,14), (0,15), (0,16), (0,17), (0,18)],
        label=["l1_1" "l1_2" "l1_3" "l2_1" "l2_2" "l2_3" "l3_1" "l3_2" "l3_3"],
        legend=:outerright,
        #c=param.color,
        #linestyle=[:solid :dash :dot],
        ylabel="l_dot  [m/s]"
    )



    fig_I = plot(
        fig0, fig1, fig2,
        layout=(3,1), size=(800, 1600),
        margin = 1.2Plots.cm,
    )
    savefig(fig_I, "sol.png")


    
    println(length(sol.t))
    println("done!")
    make_animation(sol, calc_xd)
    sol
end


@time sol = exmample();