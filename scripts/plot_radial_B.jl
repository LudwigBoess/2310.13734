using PyPlot, PyPlotUtility
using Printf, DelimitedFiles



function plot_radial_B(plot_name)

    Bfield_models = [L"B_\mathrm{sim}",
        L"B_{\beta = 50}",
        L"B_{\mathcal{F} = 0.1}",
        L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}",
        L"B_\mathrm{dyn ↑}"]

    data_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/radial_B/"

    Bfield_names = ["B_sim", "B_beta", "B_vturb", "B_FF", "B_dyn_l", "B_dyn_h"]
    folders = ["box", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_inj"]

    sim_names = ["SLOW-CR3072", "Coma", "Coma-" * L"D_\mathrm{pp}" * "-Low",
        "Coma-" * L"D_\mathrm{pp}" * "-High"]

    bonafede = readdlm("/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/bonafede_sorted.dat")

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=4.5))
    sm2.set_array([])

    alpha_ref = 1.0
    lw = 5.0
    x_pixels = 600
    legend_font_size = 30
    fig = get_figure(5.5; x_pixels)
    plot_styling!(x_pixels, axis_label_font_size=30; legend_font_size)
    gs = plt.GridSpec(1, 6, figure=fig)

    for col = 1:length(Bfield_models)

        println("col $col")
        subplot(get_gs(gs, 0, col - 1))
        ax = gca()

        axis_ticks_styling!(ax, tick_label_size=20)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim([10.0, 2.0e4])
        ax.set_ylim([1.e-1, 20.0])
        ax.set_facecolor("white")

        #if col == 3
            xlabel("Radius  " * L"r" * " [kpc]",
                fontsize=legend_font_size)
            #ax.xaxis.set_label_coords(0.8, -0.05)
        #end
        if col == 1
            ylabel("Magnetic Field  " * L"B" * " [" * L"\mu" * "G]", fontsize=legend_font_size)
            #ax.yaxis.set_label_coords(-0.2, 0.05)
        else
            ax.set_yticklabels([])
        end

        if 3 <= col <= 5
            fac = 5.0
        else
            fac = 1.0
        end

        for i_sim = 1:length(sim_names)

            fi = data_path * "/" * folders[i_sim] * "/" * Bfield_names[col] * ".dat"
            d = readdlm(fi)
            plot(d[:,1], fac .* 1.e6 .* d[:,2], 
                alpha=alpha_ref,
                color=sm2.to_rgba(i_sim), lw=lw,
                label=sim_names[i_sim])

        end

        plot(bonafede[:, 1], bonafede[:, 2], color="k", lw=lw, linestyle="--",
            label="Bonafede et al. (2010)")

        ax.text(20.0, 10.0, Bfield_models[col], fontsize=legend_font_size,
            horizontalalignment="left", verticalalignment="center")

    end

    legend(frameon=false, bbox_to_anchor=(-2.0, -0.4), ncol=5, loc="lower center")

    subplots_adjust(hspace=0.0, wspace=0.0)
    savefig(plot_name, bbox_inches="tight", transparent=false)
    close(fig)
end

map_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/coma/"
folders = ["box", "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_HB07"]
plot_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/"
plot_name = plot_path * "radial_B.pdf"

plot_radial_B(plot_name)


filenames = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]

snaps = ["036", "012", "074", "012", "074"]

Ncols = length(filenames)
Nrows = 5

α_files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_slope_$filename.xz.fits"
           for i ∈ 1:Nrows, filename ∈ filenames]

jν_files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_F_beam_1'_$filename.xz.fits"
            for i ∈ 1:Nrows, filename ∈ filenames]




folders = ["box", "zoom_dpp_1e-17", "zoom_dpp_5e-17", 
    "zoom_inj"]
data_paths = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/radial_B/" .* folders .* "/"

d = readdlm(data_paths[1] * "B_sim.dat")

bonafede = readdlm("/gpfs/work/pn68va/di67meg/Paper/LocalUniverse/synchrotron_web/data/coma_B_radial.txt")

sorted = sortperm(bonafede[:, 1])

dummy = bonafede[sorted,:]

writedlm("/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/bonafede_sorted.dat", dummy)

