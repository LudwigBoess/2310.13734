include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "shared.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using Statistics
using ProgressMeter
using SpectralCRsUtility
using Base.Threads
using PyCall
# needs to by imported by hand to make inset axis
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
@info "done"

function plot_phase_map(mass_map, x_lim, y_lim, c_lim, cmap,
    contour_map, contour_limits, contour_label,
    plot_name, transparent=false)

    if transparent
        color = "w"
    else
        color = "k"
    end

    xlabel_text = "Electron Density  " * L"n_e" * "  [cm" * L"^{-3}]"
    ylabel_text = "Temperature  " * L"T" * "  [K]"

    X = 10.0 .^ LinRange(log10(x_lim[1]), log10(x_lim[2]), size(mass_map, 1))
    Y = 10.0 .^ LinRange(log10(y_lim[1]), log10(y_lim[2]), size(mass_map, 2))

    mass_cmap = "bone_r"
    mass_label = "Mass " * L"M \:\: [M_\odot]"

    height_ratios = [0.05, 1]


    fig = get_figure(1.0, x_pixels=900)
    plot_styling!(900; color)

    gs = plt.GridSpec(2, 2, figure=fig, width_ratios=[1, 0.05]; height_ratios)

    sm = plt.cm.ScalarMappable(norm=matplotlib.colors.LogNorm(
            vmin=c_lim[1], vmax=c_lim[2]),
        cmap=mass_cmap)
    sm.set_array([])

    subplot(get_gs(gs, 0, 0))
    cax = gca()
    cb = colorbar(sm, cax=cax, orientation="horizontal")

    cb.set_label(mass_label)
    cb_ticks_styling!(cb)
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cb.ax.xaxis.set_label_coords(0.5, 2.8)

    subplot(get_gs(gs, 1, 0))
    ax = gca()

    axis_ticks_styling!(ax; color)
    ax.set_xscale("log")
    ax.set_yscale("log")

    ylabel(ylabel_text)
    xlabel(xlabel_text)

    cmap = plt.get_cmap(cmap)
    cmap.set_bad("white")

    # mass
    im1 = pcolormesh(X, Y, mass_map,
        cmap=mass_cmap,
        norm=matplotlib.colors.LogNorm(vmin=c_lim[1], vmax=c_lim[2]))

    cont1 = contour(X, Y, contour_map,
        levels=10.0 .^ LinRange(log10(contour_limits[1]),
            log10(contour_limits[2]), 20),
        cmap=cmap,
        norm=matplotlib.colors.LogNorm(vmin=contour_limits[1], vmax=contour_limits[2]),
        alpha=0.5)

    # define inset axis with zoomed window
    axins = inset_locator.inset_axes(ax,
        width="100%", height="100%",
        bbox_to_anchor=(0.62, 0.035, 0.4, 0.4),
        bbox_transform=ax.transAxes)


    axins.set_xticklabels([])
    axins.set_yticklabels([])
    axins.set_xticks([])
    axins.set_yticks([])
    axins.imshow(mass_map,
        norm=matplotlib.colors.LogNorm(vmin=c_lim[1], vmax=c_lim[2]),
        cmap=mass_cmap,
        origin="lower"
    )


    subplot(get_gs(gs, 1, 1))

    cax = gca()
    # cb = colorbar(cont1, cax=cax, use_gridspec=true)
    sm = plt.cm.ScalarMappable(norm=matplotlib.colors.LogNorm(
            vmin=contour_limits[1], vmax=contour_limits[2]),
        cmap=cmap)
    sm.set_array([])
    cb = colorbar(sm, cax=cax, use_gridspec=true)
    cb.set_label(contour_label)
    cb_ticks_styling!(cb)
    #cb.ax.yaxis.set_label_coords(3.0, 0.5)


    subplots_adjust(hspace=0.05, wspace=0.05)

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end



output_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/phase_maps/"
filename = output_path * "phase_map_mass.dat"
x_lim, y_lim, phase_map = read_phase_map(filename)

filename = output_path * "phase_map_CReE_high.dat"
x_lim, y_lim, contour_map = read_phase_map(filename)

c_lim = [1.e11, 1.e16]
contour_limits = [1.e53, 1.e59]
contour_label = L"E_\mathrm{CR,e > 1 GeV}" * " [erg]"

plot_name = plot_path * "Fig10.pdf"

plot_phase_map(phase_map, x_lim, y_lim, c_lim, "Blues_r", contour_map, 
                contour_limits, contour_label, plot_name)
