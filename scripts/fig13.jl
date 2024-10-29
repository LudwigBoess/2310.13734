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

function plot_shob_phase_map_compare(fi1, fi2, c_lim,
    plot_name; transparent=false)

    if transparent
        color = "w"
    else
        color = "k"
    end

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=4.5))
    sm2.set_array([])

    x_lim, y_lim, mass_map = read_phase_map(fi1)
    #mass_map = phase_map .* 8.479961392950551e7

    xlabel_text = "Electron Density  " * L"n_e" * "  [cm" * L"^{-3}]"
    ylabel_text = "Shock Obliquity  " * L"\theta_B" * "  [" * L"^\circ" * "]"

    X = 10.0 .^ LinRange(x_lim[1], x_lim[2], size(mass_map, 1))
    Y = LinRange(y_lim[1], y_lim[2], size(mass_map, 2))

    mass_cmap = "bone_r"
    mass_label = "Mass " * L"M \:\: [M_\odot]"

    height_ratios = [0.05, 1]


    fig = get_figure(1.0, x_pixels=900)
    plot_styling!(900; color)

    gs = plt.GridSpec(2, 1, figure=fig; height_ratios)

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
    #ax.set_xlim([X[1], X[end]])
    ax.set_xlim([2.e-9, 2.e-1])
    ax.set_ylim([Y[1], Y[end]])

    ylabel(ylabel_text)
    xlabel(xlabel_text)

    # mass
    im1 = pcolormesh(X, Y, mass_map,
        cmap=mass_cmap,
        norm=matplotlib.colors.LogNorm(vmin=c_lim[1], vmax=c_lim[2]))

    # get peak values 
    x_lim, y_lim, phase_map = read_phase_map(fi2)

    θ_mean = Vector{Float64}(undef, size(phase_map, 2))
    for i = 1:length(θ_mean)
        θ_mean[i] = Y[findmax(phase_map[:, i])[2]]
    end

    θ_mean[X.<1.e-8] .= NaN
    θ_mean[X.>1.e-1] .= NaN
    plot(X, θ_mean, color=sm2.to_rgba(4), lw=4, label="Coma")

    # get peak values 
    θ_mean = Vector{Float64}(undef, size(mass_map, 2))
    for i = 1:length(θ_mean)
        θ_mean[i] = Y[findmax(mass_map[:, i])[2]]
    end

    θ_mean[X.<1.e-8] .= NaN
    θ_mean[X.>1.e-1] .= NaN
    plot(X, θ_mean, color=sm2.to_rgba(2), lw=4, label="SLOW-CR3072" * L"^3")

    legend(frameon=false)

    subplots_adjust(hspace=0.05, wspace=0.05)

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end

data_path = "/Users/ludwigboess/Documents/Code/PaperRepos/2310.13734/data/"
c_lim = [1.e9, 3.e12]
fi1 = data_path * "phase_maps/box/phase_map_shob.dat"
fi2 = data_path * "phase_maps/zoom_inj/phase_map_shob.dat"
x_lim, y_lim, mass_map = read_phase_map(fi1)
#mass_map = phase_map .* 8.479961392950551e7

maximum(mass_map)

plot_name = plot_path * "Fig13.pdf"

plot_shob_phase_map_compare(fi1, fi2, c_lim,
    plot_name, transparent=false)
