include(joinpath(@__DIR__, "config.jl"))

using Printf
using DelimitedFiles
using PyPlot, PyPlotUtility


function read_histograms(data_file)
    d = readdlm(data_file)
    return d[:, 1], d[:, 2]./sum(d[:, 2])
end

function plot_histograms(sim_names, Bfield_names, plot_name)

    Bfield_models = [L"B_\mathrm{sim}",
                    L"B_{\beta = 50}",
                    L"B_{\mathcal{F}}",
                    L"B_\mathrm{ff}",
                    L"B_\mathrm{dyn ↓}",
                    L"B_\mathrm{dyn ↑}"]

    sim_labels = ["SLOW-CR3072" * L"^3",
                    "Coma",
                    "Coma-" * L"D_{\mathrm{pp}}"]

    Nrows = length(sim_names)

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=6.5))
    sm2.set_array([])

    x_pixels = 1200
    fig = get_figure(0.33; x_pixels)
    plot_styling!(x_pixels, axis_label_font_size=6)
    gs = plt.GridSpec(Nrows, 1, figure=fig)

    for i = 0:Nrows-1

        ax = subplot(get_gs(gs, i, 0))
        axis_ticks_styling!(ax)
        ax.set_yscale("log")
        ax.set_xlim(-12.0, -4.0)
        ax.set_ylim(1e-4, 3.0)

        if i == Nrows-1
            ax.set_xlabel(L"\log_{10}(B)" * " [" * L"G" * "]")
        else
            ax.set_xticklabels([])
        end
        ax.set_ylabel(L"N/N_\mathrm{tot}")

        axvline(log10(3.2e-6), color="k", linestyle="--"#, label=L"B_\mathrm{CMB}"
            )
        text(-11.5, 1.0, sim_labels[i+1], fontsize=16, horizontalalignment="left")

        for j = 1:length(Bfield_models)
            sim_name = sim_names[i+1] * "_" * Bfield_names[j]
            data_file = data_path * "B_histograms/$(sim_name)_B_histograms.txt"
            data = readdlm(data_file)
            plot(data[:,1], data[:,2] ./ sum(data[:,2]), label=Bfield_models[j],
                color=sm2.to_rgba(j), lw=2)
        end

        ax.xaxis.set_major_locator(plt.LinearLocator(5))


        if i == 0
            ax.legend(fontsize=16, frameon=false, bbox_to_anchor=(0.5, 1.025), ncol=3, loc="lower center")
        end

    end

    subplots_adjust(hspace=0.0, wspace=0.0)
    savefig(plot_name, bbox_inches="tight", transparent=false)
    close(fig)
end


sim_names = ["box", "zoom_inj", "zoom_dpp"]
Bfield_names = ["sim", "ff", "beta", "vturb1", "dyn_l", "dyn_h"]
plot_name = "/Users/ludwigboess/Documents/Code/PaperRepos/2310.13734/Plots/Fig14.pdf"

plot_histograms(sim_names, Bfield_names, plot_name)