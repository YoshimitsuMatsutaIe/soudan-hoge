"""シミュレーション"""


using LinearAlgebra  # 逆行列とノルムの計算に使用
using Plots  # グラフ作成
#using ProgressBars
#using Parameters
using DifferentialEquations  # 微分方程式のソルバ
using LaTeXStrings


include("kinematics.jl")
include("differential_kinematics.jl")
#include("controller.jl")
include("plot_utils.jl")

using .Kinematics
#using .Controller



const R = 0.1

"""適当な目標位置変位"""
function calc_xd(t::T) where T
    [
        R * sin(3t),
        R * cos(3t),
        0.147,
    ]
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



"""ソルバで使うやつ"""
function state_eq!(
    q_dot::Vector{T}, q::Vector{T},
    p,
    t::T
    ) where T
    q = X[1:3]
    q_dot = X[4:6]
    H = X[7:9]
    #println(t)
    τ = calc_torque(
    p, q, q_dot,
    calc_qd(t), calc_qd_dot(t), calc_qd_dot_dot(t),
    )
    X_dot[1:3] = q_dot
    X_dot[4:6] = calc_q_dot_dot(τ, q, q_dot, H)
    X_dot[7:9] = H_dot(H, q, q_dot)
end


"""微分方程式を解く"""
function run_simulation(
    TIME_SPAN::T, method_name::String, controller_param
    ) where T

    X₀ = zeros(T, 9)
    t_span = (0.0, TIME_SPAN) 
    println(method_name * " now...")
    prob = ODEProblem(state_eq!, X₀, t_span, controller_param)
    solve(prob)
end



function reproduce_τ(p::PassivityBasedAdaptiveController{T}, t::T, u::Vector{T}) where T
    calc_torque(
        p, u[1:3], u[4:6],
        calc_qd(t), calc_qd_dot(t),
        u[10:15]
    )
end


function reproduce_τ(p, t::T, u::Vector{T}) where T
    calc_torque(
        p, u[1:3], u[4:6],
        calc_qd(t), calc_qd_dot(t), calc_qd_dot_dot(t),
    )
end


"""実行"""
function exmample()

    TIME_SPAN = 0.05

    hutashikasa = false  # 剛性行列と減衰行列に不確かさがあるかないか

    params = [
        # (
        #     name = "kine",
        #     p = KinematicController(
        #         K_kin = Dynamics.K,
        #         isUncertainty = hutashikasa,
        #     ),
        #     color = :red
        # ),
        # (
        #     name = "pdfb",
        #     p = PDandFBController(
        #         Kd = Matrix{Float64}(I, 3, 3)*200,
        #         Kp = Matrix{Float64}(I, 3, 3)*10000,
        #         isUncertainty = hutashikasa,
        #     ),
        #     color = :blue,
        # ),
        # (
        #     name = "psiv",
        #     p = PassivityBasedController(
        #         Λ = Matrix{Float64}(I, 3, 3)*10,
        #         KG = Matrix{Float64}(I, 3, 3)*10,
        #         isUncertainty = hutashikasa,
        #     ),
        #     color = :magenta
        # ),
        # (
        #     name = "psad",
        #     p = PassivityBasedAdaptiveController(
        #         invΓ = Matrix{Float64}(I, 6, 6)*1e+5,
        #         Λ = Matrix{Float64}(I, 3, 3)*10,
        #         KG = Matrix{Float64}(I, 3, 3)*10,
        #         isUncertainty = hutashikasa,
        #     ),
        #     color = :black
        # ),
        (
            name = "sdre",
            p = SDREController(
                Q = Matrix{Float64}(I, 6, 6)*2000,
                R = Matrix{Float64}(I, 3, 3)*0.1,
                isUncertainty = hutashikasa,
            ),
            color = :green
        ),
        # (
        #     name = "dmpc",
        #     p = MPCController(
        #         Q = diagm([1, 1, 1, 1, 1, 1, 0.0, 0.0, 0.0]),
        #         R = diagm([1.0, 1.0, 1.0,]),
        #         n = 10,
        #         Δt = 1.0e-15,
        #         isUncertainty = hutashikasa
        #     ),
        #     color = :cyan
        # )
    ]
    sols = []
    fig0 = plot(xlims=(0.0, TIME_SPAN), ylims=(0,0.02))
    fig_tau = plot(xlims=(0.0, TIME_SPAN), ylims=(-150.0, 150.0))
    fig1 = plot(xlims=(0.0, TIME_SPAN))
    fig2 = plot(xlims=(0.0, TIME_SPAN))
    fig3 = plot(xlims=(0.0, TIME_SPAN))
    

    for param in params
        @time sol = run_simulation(TIME_SPAN, param.name, param.p)
        push!(sols, sol)

        error = Vector{Vector{Float64}}(undef, length(sol.t))
        τ = Vector{Vector{Float64}}(undef, length(sol.t))
        for (i, t) in enumerate(sol.t)
            error[i] = calc_qd(t) .- sol.u[i][1:3]
            τ[i] = reproduce_τ(param.p, t, sol.u[i])
        end
        L2_error = norm.(error)
        plot!(
            fig0,
            sol.t, L2_error,
            label=(param.name * "-") * "e",
            legend=:outerright,
            c=param.color,
            ylabel="error [m]"
        )

        x, y, z = split_vec_of_arrays(τ)
        plot!(fig_tau, sol.t, x, label=param.name*"-"*"1", c=param.color, linestyl=:solid)
        plot!(fig_tau, sol.t, y, label=param.name*"-"*"2", c=param.color, linestyle=:dash)
        plot!(fig_tau, sol.t, z, label=param.name*"-"*"3", c=param.color, linestyle=:dot)
        plot!(
            fig_tau, legend=:outerright,
            ylabel="τ [N/rad]"
        )

        plot!(
            fig1,
            sol, vars=[(0,1), (0,2), (0,3)],
            label=(param.name * "-") .* ["1" "2" "3"],
            legend=:outerright,
            c=param.color,
            linestyle=[:solid :dash :dot],
            ylabel="l [m]"
        )
        plot!(
            fig2,
            sol, vars=[(0,4), (0,5), (0,6)],
            label=(param.name * "-") .* ["1" "2" "3"],
            legend=:outerright,
            c=param.color,
            linestyle=[:solid :dash :dot],
            ylabel="l_dot  [m/s]"
        )
        plot!(
            fig3,
            sol, vars=[(0,7), (0,8), (0,9)],
            label=(param.name * "-") .* ["1" "2" "3"],
            legend=:outerright,
            c=param.color,
            linestyle=[:solid :dash :dot],
            ylabel="h  [N/rad]"
        )
    
    end

    fig_I = plot(
        fig0, fig_tau, fig1, fig2, fig3,
        layout=(5,1), size=(800, 1600),
        margin = 1.2Plots.cm,
    )
    savefig(fig_I, "solve.png")


    sols
    println("done!")
end


@time _ = exmample()