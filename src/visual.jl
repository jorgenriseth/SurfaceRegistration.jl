using PyPlot

function plot_parametric_surface(func; n_points=41, axis_title="", axis=nothing, camera=(15, 45))
    # Plot Data
    xs = ys = range(0, 1, length=n_points)

    F = [func([xi, yi]) for xi in xs, yi in xs]
    fx = [el[1] for el in F]
    fy = [el[2] for el in F]
    fz = [el[3] for el in F]

    if axis === nothing
        fig = figure()
        plot_surface(fx, fy, fz)
        xlabel("x")
        ylabel("y")
        zlabel("z")
        title(axis_title)
        gca().view_init(camera...)
        return fig
    end

    axis.plot_surface(fx, fy, fz)
    axis.set_ylabel("y")
    axis.set_zlabel("z")
    axis.set_xlabel("x")
    axis.set_title(axis_title)
    axis.view_init(camera[1], camera[2])
    return axis
end

function plot_parametric_wireframe(func; n_points=11, axis_title="", axis=nothing, camera=(15, 45))
    # Plot Data
    xs = ys = range(0, 1, length=n_points)

    F = [func([xi, yi]) for xi in xs, yi in xs]
    fx = [el[1] for el in F]
    fy = [el[2] for el in F]
    fz = [el[3] for el in F]

    if axis === nothing
        fig = figure()
        plot_wireframe(fx, fy, fz)
        xlabel("x")
        ylabel("y")
        zlabel("z")
        title(axis_title)
        gca().view_init(camera...)
        return fig
    end

    axis.plot_wireframe(fx, fy, fz)
    axis.set_ylabel("y")
    axis.set_zlabel("z")
    axis.set_xlabel("x")
    axis.set_title(axis_title)
    axis.view_init(camera[1], camera[2])
    return axis
end

function plot_diffeomorphism(ψ; n_points=15, axis=nothing)
    xs = 0.:(1/n_points):1.
    
    if axis == nothing
        fig = figure()
        
        # Plot horizontal lines
        for xi in xs
            line_x, line_y = [xi for yi in xs], [yi for yi in xs]
            warped_line_x = [ψ([xi, yi])[1] for yi in xs]
            warped_line_y = [ψ([xi, yi])[2] for yi in xs]

            plot(line_x, line_y, c="r", lw=0.3)
            plot(warped_line_x, warped_line_y, c="k")
        end
        
        # Plot Warped Lines Vertically
        for yi in xs
            line_x, line_y = [xi for xi in xs], [yi for xi in xs]
            warped_line_x = [ψ([xi, yi])[1] for xi in xs]
            warped_line_y = [ψ([xi, yi])[2] for xi in xs]

            plot(line_x, line_y, c="r", lw=0.3)
            plot(warped_line_x, warped_line_y, c="k")
        end
        
        return fig
    end
    
    
        # Plot horizontal lines
    for xi in xs
        line_x, line_y = [xi for yi in xs], [yi for yi in xs]
        warped_line_x = [ψ([xi, yi])[1] for yi in xs]
        warped_line_y = [ψ([xi, yi])[2] for yi in xs]

        axis.plot(line_x, line_y, c="r", lw=0.3)
        axis.plot(warped_line_x, warped_line_y, c="k")
    end

    # Plot Warped Lines Vertically
    for yi in xs
        line_x, line_y = [xi for xi in xs], [yi for xi in xs]
        warped_line_x = [ψ([xi, yi])[1] for xi in xs]
        warped_line_y = [ψ([xi, yi])[2] for xi in xs]

        axis.plot(line_x, line_y, c="r", lw=0.3)
        axis.plot(warped_line_x, warped_line_y, c="k")
    end
end