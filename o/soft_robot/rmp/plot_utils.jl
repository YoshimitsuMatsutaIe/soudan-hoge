"""視覚化"""

using Plots
using ProgressBars


include("utils.jl")
include("kinematics.jl")

using .Kinematics


function make_plot_basic(data)
    x, y, z = split_vec_of_arrays(data.q)
    fig = plot(xlim=(0, data.t[end]))
    plot!(fig, data.t, x, label="l1")
    plot!(fig, data.t, y, label="l2")
    plot!(fig, data.t, z, label="l3")

    x, y, z = split_vec_of_arrays(data.q_dot)
    fig2 = plot(xlim=(0, data.t[end]))
    plot!(fig2, data.t, x, label="l1_dot")
    plot!(fig2, data.t, y, label="l2_dot")
    plot!(fig2, data.t, z, label="l3_dot")

    x, y, z = split_vec_of_arrays(data.error)
    fig3 = plot(xlim=(0, data.t[end]))
    plot!(fig3, data.t, x, label="e1")
    plot!(fig3, data.t, y, label="e2")
    plot!(fig3, data.t, z, label="e3")

    x, y, z = split_vec_of_arrays(data.τ)
    fig4 = plot(xlim=(0, data.t[end]))
    plot!(fig4, data.t, x, label="τ1")
    plot!(fig4, data.t, y, label="τ2")
    plot!(fig4, data.t, z, label="τ3")

    fig_I = plot(
        fig, fig2, fig3, fig4,
        layout=(4, 1), size=(500, 800))
    savefig(fig_I, "julia.png")

end







"""1フレームを描写"""
function draw_frame(
    t::T, q::Vector{T}, qd::Union{Vector{T}, Nothing},
    fig_shape
    ) where T
    

    arm = Arm(q)

    x, y, z = split_vec_of_arrays(arm)
    fig = plot(
        x, y, z,
        #marker=:circle,
        aspect_ratio = 1,
        #markersize=2,
        label="arm",
        xlabel = "X[m]", ylabel = "Y[m]", zlabel = "Z[m]",
    )

    if !isnothing(qd)
        scatter!(
            fig,
            [xd[1]], [xd[2]], [xd[3]],
            label="qd",
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
        # xlims=(fig_shape.xl, fig_shape.xu),
        # ylims=(fig_shape.yl, fig_shape.yu),
        # zlims=(fig_shape.zl, fig_shape.zu),
        legend = true,
        size=(600, 600),
        title = string(round(t, digits=2)) * "[s]"
    )

    return fig, arm
end


fig, arm = draw_frame(
    0.0,
    [
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.05,
        0.0, 0.05, 0.0
    ], nothing, nothing)

"""アニメ作成"""
function make_animation(data)
    println("アニメ作成中...")
    # 枚数決める
    #println(data.t)
    epoch_max = 100
    epoch = length(data.t)
    if epoch < epoch_max
        step = 1
    else
        step = div(epoch, epoch_max)
    end

    #println(step)

    x_max = 0.1
    x_min = -0.1
    y_max = 0.1
    y_min = -0.1
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

    xd_all = [Phi0(q, 1.0) for q in data.qd]
    x, y, z = split_vec_of_arrays(xd_all)

    anim = Animation()
    @gif for i in tqdm(1:step:length(data.t))
        _fig = draw_frame(data.t[i], data.q[i], data.qd[i], fig_shape)
        plot!(_fig, x, y, z, label="xd",)
        frame(anim, _fig)
    end



    gif(anim, "julia.gif", fps = 60)

    println("アニメ作成完了")
end
