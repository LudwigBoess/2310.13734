include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "shared.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Statistics
using Printf
using ProgressMeter
using PyCall
# needs to by imported by hand to make inset axis
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
axes_divider = pyimport("mpl_toolkits.axes_grid1.axes_divider")
@info "done"

function plot_phase_maps(phase_map_path, plot_name)

    Bfield_models = [L"B_\mathrm{sim}",
        L"B_{\beta}",
        L"B_{\mathcal{F}}",
        L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}",
        L"B_\mathrm{dyn ↑}"]

    filename = phase_map_path .* ["box/phase_map_synch_emissivity_144MHz_$B.dat"
                             for B ∈ ["B_sim", "B_beta50", "B_01Pturb",
        "B_FF", "B_dyn_l", "B_dyn_h"]
    ]

    filename_1D = phase_map_path .* ["box/bin_1D_synch_emissivity_144MHz_$B.dat"
                                for B ∈ ["B_sim", "B_beta50", "B_01Pturb",
        "B_FF", "B_dyn_l", "B_dyn_h"]
    ]

    filename_1D_zoom = phase_map_path .* ["zoom_inj/bin_1D_synch_emissivity_144MHz_$B.dat"
                                     for B ∈ ["B_sim", "B_beta50", "B_01Pturb",
        "B_FF", "B_dyn_l", "B_dyn_h"]
    ]

    filename_1D_zoom_Dpp_low = phase_map_path .* ["zoom_dpp_low/bin_1D_synch_emissivity_144MHz_$B.dat"
                                             for B ∈ ["B_sim", "B_beta50", "B_01Pturb",
        "B_FF", "B_dyn_l", "B_dyn_h"]
    ]


    filename_1D_zoom_Dpp_high = phase_map_path .* ["zoom_dpp_high/bin_1D_synch_emissivity_144MHz_$B.dat"
                                              for B ∈ ["B_sim", "B_beta50", "B_01Pturb",
        "B_FF", "B_dyn_l", "B_dyn_h"]
    ]


    x_lim, y_lim, phase_map_jν = read_phase_map(filename[1])

    lw = 3
    xlabel_text = "Electron Density  " * L"n_e" * "  [cm" * L"^{-3}]"
    ylabel_text = "Synch. Emissivity  " * L"j_\nu" * "  [erg s" * L"^{-1}" * " Hz" * L"^{-1}" * " cm" * L"^{-3}]"

    X = 10.0 .^ LinRange(log10(x_lim[1]), log10(x_lim[2]), size(phase_map_jν, 1))
    Y = 10.0 .^ LinRange(log10(y_lim[1]), log10(y_lim[2]), size(phase_map_jν, 2))

    mass_cmap = "bone_r"
    mass_label = "Mass " * L"M \:\: [M_\odot]"
    c_lim = [1.e7, 1.e14]

    color = "k"
    alpha_ref = 1.0

    height_ratios = [0.05, 1, 1]

    x_pixels = 1_000
    legend_font_size = 12
    fig = get_figure(1.5; x_pixels)
    plot_styling!(x_pixels; color, legend_font_size)

    gs = plt.GridSpec(3, 3, figure=fig; height_ratios)

    sm = plt.cm.ScalarMappable(norm=matplotlib.colors.LogNorm(
            vmin=c_lim[1], vmax=c_lim[2]),
        cmap=mass_cmap)
    sm.set_array([])

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=4.5))
    sm2.set_array([])

    subplot(get_gs(gs, 0, 0:3))
    cax = gca()
    cb = colorbar(sm, cax=cax, orientation="horizontal")

    cb.set_label(mass_label)
    cb_ticks_styling!(cb)
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cb.ax.xaxis.set_label_coords(0.5, 2.8)


    Nfile = 1

    for row = 1:2, col = 0:2

        println("row $row, col $col")
        subplot(get_gs(gs, row, col))
        ax = gca()

        axis_ticks_styling!(ax; color)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim([1.e-9, 1.0])
        ax.set_ylim([1.e-52, 1.e-37])
        ax.set_facecolor("white")

        if row == 1
            ax.set_xticklabels([])
        else
            if col == 1
                xlabel(xlabel_text)
            end
        end

        if col == 0 && row == 1
            ylabel(ylabel_text)
            ax.yaxis.set_label_coords(-0.2, 0.05)
        end

        cmap = plt.get_cmap(mass_cmap)
        cmap.set_bad("white")

        ax.grid(true)

        # mass
        x_lim, y_lim, phase_map_jν = read_phase_map(filename[Nfile])
        phase_map_jν .*= 8.479961392950551e7 # multiply with particle mass 
        im1 = pcolormesh(X, Y, phase_map_jν,
            cmap=mass_cmap,
            norm=matplotlib.colors.LogNorm(vmin=c_lim[1], vmax=c_lim[2]))

        if Nfile > 0
            ne_bins, jnu_mean = read_1D_data(filename_1D[Nfile])
            sel = findall(ne_bins .< 5.e-9)
            jnu_mean[sel] .= NaN
            sel = findall(ne_bins .> 3.e-1)
            jnu_mean[sel] .= NaN
            plot(ne_bins, jnu_mean, alpha=alpha_ref,
                color=sm2.to_rgba(1), lw=lw, label="SLOW-CR3072" * L"^3")

            ne_bins, jnu_mean = read_1D_data(filename_1D_zoom[Nfile])
            sel = findall(ne_bins .< 5.e-9)
            jnu_mean[sel] .= NaN
            sel = findall(ne_bins .> 3.e-1)
            jnu_mean[sel] .= NaN
            plot(ne_bins, jnu_mean, alpha=alpha_ref,
                color=sm2.to_rgba(2), lw=lw, label="Coma")#L"\textsc{Coma}")

            ne_bins, jnu_mean = read_1D_data(filename_1D_zoom_Dpp_low[Nfile])
            sel = findall(ne_bins .< 5.e-9)
            jnu_mean[sel] .= NaN
            sel = findall(ne_bins .> 3.e-1)
            jnu_mean[sel] .= NaN
            plot(ne_bins, jnu_mean, alpha=alpha_ref,
                color=sm2.to_rgba(3), lw=lw,
                label="Coma-" * L"D_\mathrm{pp}" * "-Low")

            ne_bins, jnu_mean = read_1D_data(filename_1D_zoom_Dpp_high[Nfile])
            sel = findall(ne_bins .< 5.e-9)
            jnu_mean[sel] .= NaN
            sel = findall(ne_bins .> 3.e-1)
            jnu_mean[sel] .= NaN
            plot(ne_bins, jnu_mean, alpha=alpha_ref,
                color=sm2.to_rgba(4), lw=lw,
                label="Coma-" * L"D_\mathrm{pp}" * "-High")
        end


        text(1.e-8, 1.e-39, Bfield_models[Nfile], color="k")

        errorbar([1.e-5], [1.e-44], xerr=[5.e-6], yerr=[7.e-45], uplims=true, ecolor="k", color="k")

        Nfile += 1
    end

    subplot(get_gs(gs, 1, 1))
    ax = gca()
    legend(frameon=false, bbox_to_anchor=(0.5, -1.5), ncol=4, loc="lower center")

    subplots_adjust(hspace=0.05)
    savefig(plot_name, bbox_inches="tight", transparent=false)
    close(fig)
end

phase_map_path = data_path * "phase_maps_new/"
plot_name = plot_path * "Fig08.pdf"

plot_phase_maps(phase_map_path, plot_name)

