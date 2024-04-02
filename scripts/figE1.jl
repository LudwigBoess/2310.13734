using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using ColorSchemes



function plot_it()

# found by selecting particles by hand
    coma_center = [244977.58, 327853.62, 245989.75] .* 1.e-3
    cylinder = GadgetCylinder([247.75, 332.0, 248.5] - coma_center, 
                              [248.5, 332.75, 249.25] - coma_center, 
                              2.5e3)

    files = [map_path * "coma/zoom_dpp_5e-17/coma_20Mpc_012.synch_F_beam_1'_144MHz_dyn_h.xz.fits"]

    x_pixels = 800
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

        #ax.plot([cylinder.pos_start[2], cylinder.pos_end[2]], [cylinder.pos_start[3], cylinder.pos_end[3]], c="red", alpha=1.0)

    subplots_adjust(hspace=0.0, wspace=0.12)
    plot_name = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/FigE1.pdf"

    println("saving")
    savefig(plot_name, bbox_inches="tight", transparent=false, dpi=400)

    close(fig)
end

plot_it()