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
using .DifferentialKinematics
#using .Controller



const R = 0.1
const offset_z = 0.5

"""適当な目標位置変位"""
function calc_xd(t::T) where T
    [
        R * sin(3t),
        R * cos(3t),
        offset_z,
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


function state_eq!(
    X_dot::Vector{T}, X::Vector{T},
    p,
    t::T
    ) where T

    J = J_2(X[1:9], 1.0)
    X_dot[1:9] = X[10:18]
    X_dot[10:18] = pinv(J) *
    (
        calc_xd_dot_dot(t) .-
        p.D*(J*X[10:18] .- calc_xd_dot(t)) .-
        p.K*(Phi_2(X[1:9], 1.0) - calc_xd(t))
    )
end


"""微分方程式を解く"""
function run_simulation(
    TIME_SPAN::T
    ) where T

    X₀ = zeros(T, 18)
    t_span = (0.0, TIME_SPAN)
    p = (
        K = Matrix{T}(I, 3, 3)*10,
        D = Matrix{T}(I, 3, 3)*1,
    )
    prob = ODEProblem(state_eq!, X₀, t_span, p)
    solve(prob)
end


"""実行"""
function exmample()

    TIME_SPAN = 3.0

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
        ylabel="error [m]"
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


    sol
    println("done!")
end


@time _ = exmample()