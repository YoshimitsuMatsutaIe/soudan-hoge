"""視覚化"""

using Plots


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
