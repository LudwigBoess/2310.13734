using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf


function plot_corners(corner_lower_left, corner_upper_right, plot_name)

    files = ["/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/coma/zoom_dpp_5e-17/coma_20Mpc_012.synch_F_beam_1'_dyn_h.xz.fits"]

    vmin = 1.e-8
    vmax = 1.e2

    # vmin = 1.e48
    # vmax = 1.e51
    fig = get_figure(1.0, x_pixels=800)
    plot_styling!()
    ax = gca()
    axis_ticks_styling!(ax, color="w")
    map, par = PyPlotUtility.read_map_par(1, 1, files, nothing, nothing)
    extent = [par.x_lim[1], par.x_lim[2], par.z_lim[1], par.z_lim[2]]
    x = 1
    y = 2
    xlabel("x")
    ylabel("z")


    im = ax.imshow(map, norm=matplotlib.colors.LogNorm(; vmin, vmax),
        cmap="Purples_r",
        origin="lower",
        extent=extent
    )

    plot([corner_lower_left[x], corner_upper_right[x], corner_upper_right[x], corner_lower_left[x], corner_lower_left[x]],
        [corner_lower_left[y], corner_lower_left[y], corner_upper_right[y], corner_upper_right[y], corner_lower_left[y]],
        color="gray")

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end

GU = GadgetPhysical(read_header("/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/mhd_cr8p24eDpp_5e-17/snapdir_012/snap_012.0"))

cylinder = GadgetCylinder([247.75, 332.0, 248.5] .* 1.e3,
    [248.5, 332.75, 249.25] .* 1.e3,
    2.5e3)

corner_lower_left = [247.75, 248.5] .* 1.e3 .* GU.x_physical
corner_upper_right = [248.5, 249.25] .* 1.e3 .* GU.x_physical

plot_name = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/relic_selection.png"
plot_corners(corner_lower_left, corner_upper_right, plot_name)

cube = GadgetCube(corner_lower_left, corner_upper_right)

println(GadgetIO.get_geometry_center(cube) .* 1.e-3)