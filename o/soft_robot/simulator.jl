"""シミュレーション"""


using LinearAlgebra  # 逆行列とノルムの計算に使用
using Plots  # グラフ作成
#using ProgressBars
#using Parameters
using DifferentialEquations  # 微分方程式のソルバ
#using LaTeXStrings


include("kinematics.jl")
include("dynamics.jl")
include("controller.jl")
include("plot_utils.jl")

using .Kinematics: Phi0, Arm
using .Dynamics: H_dot, calc_q_dot_dot
using .Controller



const R = 0.01

"""適当な目標アクチュエータ変位"""
function calc_qd(t::T) where T
    [
        R * sin(3t),
        R * sin(3t + π/2),
        R * sin(3t + π),
    ]
end


"""適当な目標アクチュエータ速度"""
function calc_qd_dot(t::T) where T
    [
        R * 3cos(3t),
        R * 3cos(3t + π/2),
        R * 3cos(3t + π),
    ]
end


"""適当な目標アクチュエータ加速度"""
function calc_qd_dot_dot(t::T) where T
    [
        R * -9sin(3t),
        R * -9sin(3t + π/2),
        R * -9sin(3t + π),
    ]
end



"""ソルバで使うやつ"""
function state_eq!(X_dot::Vector{T}, X::Vector{T}, p, t::T) where T
    q = X[1:3]
    q_dot = X[4:6]
    H = X[7:9]
    #println(t)
    τ = calc_torque(
    p, q, q_dot,
    calc_qd(t), calc_qd_dot(t), calc_qd_dot_dot(t)
    )

    X_dot[1:3] = q_dot
    X_dot[4:6] = calc_q_dot_dot(τ, q, q_dot, H)
    X_dot[7:9] = H_dot(H, q, q_dot)

end


"""ソルバ使用

硬い微分方程式でも正確に解いてくれる
"""
function run_simulation(;TIME_SPAN::T=1.0, method) where T

    X₀ = zeros(T, 9)
    t_span = (0.0, TIME_SPAN)
    
    parames = (
        (
            name = "kine",
            p = KinematicController(Dynamics.K),
            marker = :solid,
        ),
        (
            name = "pdfb",
            p = PDandFBController(
                Matrix{Float64}(I, 3, 3)*200,
                Matrix{Float64}(I, 3, 3)*10000#1500
            ),
            marker = :dash
        ),
        (
            name = "psiv",
            p = PassiveController(
                Matrix{Float64}(I, 3, 3)*10,
                Matrix{Float64}(I, 3, 3)*10,
            ),
            marker = :dot
        )
    )

    sols = []
    fig0 = plot(xlims=(0.0, TIME_SPAN))
    fig_tau = plot(xlims=(0.0, TIME_SPAN), ylims=(-50.0, 150.0))
    fig1 = plot(xlims=(0.0, TIME_SPAN))
    fig2 = plot(xlims=(0.0, TIME_SPAN))
    fig3 = plot(xlims=(0.0, TIME_SPAN))
    
    for param in parames
        println(param.name * " now...")
        prob = ODEProblem(state_eq!, X₀, t_span, param.p)
        sol = solve(prob)
        push!(sols, sol)

        error = Vector{Vector{T}}(undef, length(sol.t))
        τ = Vector{Vector{T}}(undef, length(sol.t))
        for (i, t) in enumerate(sol.t)
            error[i] = calc_qd(t) .- sol.u[i][1:3]
            τ[i] = calc_torque(
                param.p, sol.u[i][1:3], sol.u[i][4:6],
                calc_qd(t), calc_qd_dot(t), calc_qd_dot_dot(t)
            )
        end
        L2_error = norm.(error)
        plot!(
            fig0,
            sol.t, L2_error,
            label=(param.name * "-") * "err",
            legend=:outerright,
            linestyle=param.marker
        )

        x, y, z = split_vec_of_arrays(τ)
        plot!(fig_tau, sol.t, x, label=param.name*"-"*"τ1_", linestyle=param.marker)
        plot!(fig_tau, sol.t, y, label=param.name*"-"*"τ2_", linestyle=param.marker)
        plot!(fig_tau, sol.t, z, label=param.name*"-"*"τ3_", linestyle=param.marker)
        plot!(fig_tau,legend=:outerright)
    
        plot!(
            fig1,
            sol, vars=[(0,1), (0,2), (0,3)],
            label=(param.name * "-") .* ["l1_" "l2_" "l3_"],
            legend=:outerright,
            linestyle=param.marker
        )
        plot!(
            fig2,
            sol, vars=[(0,4), (0,5), (0,6)],
            label=(param.name * "-") .* ["dl1" "dl2" "dl3"],
            legend=:outerright,
            linestyle=param.marker
        )
        plot!(
            fig3,
            sol, vars=[(0,7), (0,8), (0,9)],
            label=(param.name * "-") .* ["h1_" "h2_" "h3_"],
            legend=:outerright,
            linestyle=param.marker
        )

    end
    fig_I = plot(
        fig0, fig_tau, fig1, fig2, fig3,
        layout=(5,1), size=(500, 1200)
    )
    savefig(fig_I, "solve.png")
    sols
end



@time sol = run_simulation(TIME_SPAN=1.0, method=nothing)