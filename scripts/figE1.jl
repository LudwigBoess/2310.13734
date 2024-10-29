using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using ColorSchemes

include(joinpath(@__DIR__, "config.jl"))

function plot_it()

# found by selecting particles by hand
    # coma_center = [244977.58, 327853.62, 245989.75] .* 1.e-3
    # cylinder = GadgetCylinder([247.75, 332.0, 248.5] - coma_center, 
    #                           [248.5, 332.75, 249.25] - coma_center, 
    #                           2.5e3)

    files = [map_path * "zoom_inj/coma_20Mpc_074.synch_F_beam_1'_144MHz_dyn_h.xz.fits"]

    x_pixels = 900
    fig = get_figure(1.0; x_pixels)
    plot_styling!(x_pixels)

        ax = gca()
        axis_ticks_styling!(ax, color="white")
        xlabel(L"x" * " [Mpc]")
        ylabel(L"y" * " [Mpc]")

        map, par = PyPlotUtility.read_map_par(1, 1, files, nothing, nothing)

        norm = matplotlib.colors.LogNorm(vmin=1.e-8, vmax=1.e2)
            im = get_imshow(ax, copy(map),
            1.e-3 .* par.x_lim .- (par.center[1] * 1.e-3),
            1.e-3 .* par.y_lim .- (par.center[2] * 1.e-3),
            cnorm=norm, cmap=ColorMap("bukavo", ColorSchemes.bukavu[1:end]))

        cb = get_colorbar_top(ax, im, "Synchrotron Flux  " * L"F_{ν,144 \mathrm{ MHz}}" * " [mJy arcmin" * L"^{-2}" * "]")
        locmaj = matplotlib.ticker.LogLocator(base=10, numticks=12)
        cb.xaxis.set_major_locator(locmaj)
        locmin = plt.LogLocator(base=10.0, subs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks=20)
        cb.xaxis.set_minor_locator(locmin)
        cb.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        cb.xaxis.set_label_coords(0.5, 2.0 + 0.06)
        #ax.plot([cylinder.pos_start[2], cylinder.pos_end[2]], [cylinder.pos_start[3], cylinder.pos_end[3]], c="red", alpha=1.0)

    subplots_adjust(hspace=0.0, wspace=0.12)
    plot_name = plot_path * "FigE1.png"

    println("saving")
    savefig(plot_name, bbox_inches="tight", transparent=false, dpi=400)

    close(fig)
end

plot_it()