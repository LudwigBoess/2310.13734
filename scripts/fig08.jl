include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "shared.jl"))

using PyPlot, PyPlotUtility
using Printf
using DelimitedFiles

function plot_spectra(spectra_path, folders, sim_names, Bfield_names, plot_name)

    Bfield_models = [L"B_\mathrm{sim}",
        L"B_{\beta}",
        L"B_{\mathcal{F}}",
        L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}",
        L"B_\mathrm{dyn ↑}"]

    Nrows = 2
    Ncols = length(folders)

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=6.5))
    sm2.set_array([])

    x_pixels = 1100
    fig = get_figure(1.5; x_pixels)
    plot_styling!(x_pixels, axis_label_font_size=14)
    gs = plt.GridSpec(Nrows, Ncols, figure=fig)

    for i = 0:Ncols-1

        # energy spectrum
        ax = subplot(get_gs(gs, 0, i))
        axis_ticks_styling!(ax)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(8, 2e5)
        # ax.set_ylim(1e-3, 1e2)
        # ax.set_ylim(1e9, 1e14)
        ax.set_ylim(1e-51, 1e-30)

        ax.set_xlabel(L"\hat{p}" * " [" * L"(m_e c)^{-1}" * "]")

        if i == 0
            # ax.set_ylabel(L"N_e" * " [arb. units]")
            ax.set_ylabel(L"f(p)" * " [arb. units]")
        else
            ax.set_yticklabels([])
        end

        filename = spectra_path * "$(folders[i+1])/f_spectrum.dat"
        data = readdlm(filename)
        plot(data[:, 1], data[:, 2], color="darkblue", lw=3)
        text(1.e5, 1.0e-32, sim_names[i+1], fontsize=20, horizontalalignment="right")

        plot([4.e2, 4.e3], [1.e-35, 1.e-39], color="k", lw=2, linestyle="--")
        text(7.e2, 3.e-38, L"q = -4", fontsize=20, rotation=-35)

        get_cr_energy_axis!(ax, "e")

        # synchrotron spectrum
        ax = subplot(get_gs(gs, 1, i))
        axis_ticks_styling!(ax)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(8.0, 2.e4)
        ax.set_ylim(1e15, 1e20)
        ax.set_xlabel(L"\nu" * " [MHz]")

        if i == 0
            ax.set_ylabel(L"P_\nu" * " [W Hz" * L"^{-1}" * "]")
        else
            ax.set_yticklabels([])
        end

        for j = 1:length(Bfield_models)
            filename = spectra_path * "$(folders[i+1])/synch_spectrum_$(Bfield_names[j]).dat"
            data = readdlm(filename)
            plot(data[:, 1] .* 1.e-6, data[:, 2], label=Bfield_models[j],
                color=sm2.to_rgba(j), lw=3)
        end

        P0 = 1.e19
        plot([5.e1, 5.e2], P0 .* [1.0, 10.0^(-1.5)], color="k", lw=3, linestyle="--")
        text(9.e1, 1.5e18, L"\alpha_\nu = -1.5", fontsize=20, rotation=-40)
    end

    handles, labels = gca().get_legend_handles_labels()
    #order = [1, 6, 2, 3, 4, 5]
    order = collect(1:6)
    l = legend([handles[idx] for idx in order], [Bfield_models[idx] for idx in order],
        frameon=false, bbox_to_anchor=(-0.5, -0.5), ncol=6, loc="lower center")

    subplots_adjust(hspace=0.3, wspace=0.0)
    savefig(plot_name, bbox_inches="tight", transparent=false)
    close(fig)
end

spectra_path = data_path * "spectra/"
folders = ["box", "zoom_inj", "zoom_dpp"]
Bfield_names = ["sim", "beta", "vturb1", "ff", "dyn_l", "dyn_h"]

sim_names = ["SLOW-CR3072" * L"^3",
             "Coma",
             "Coma-" * L"D_{\mathrm{pp}}"]

plot_name = plot_path * "Fig08.pdf"
plot_spectra(spectra_path, folders, sim_names, Bfield_names, plot_name)