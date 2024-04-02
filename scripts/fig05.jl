include(joinpath(@__DIR__, "config.jl"))

using PyPlot, PyPlotUtility
using Printf, DelimitedFiles

"""
    Maps
"""

using GadgetIO, GadgetUnits
using Printf 
using PyPlot, PyPlotUtility
using PyCall
cm = pyimport("cmasher")



Bfield_names = ["B_sim", "B_beta", 
    "B_vturb_rescaled",
    "B_FF_rescaled", "B_dyn_l_rescaled", "B_dyn_h_rescaled"]


files = @. map_path * "coma/zoom_dpp_1e-17/coma_5Mpc_$(@sprintf("%03i", 74)).RM_" * Bfield_names * ".a.y.fits"

Ncols = length(Bfield_names)
Nrows = 1

vmin_arr = [-300.0]
vmax_arr = [300.0]

im_cmap = ["seismic"]

cb_labels = ["Rotation Measure  " * L"\mathrm{RM}" * " [rad m" * L"^{-2}" * "]"]

annotate_time = trues(Nrows * Ncols)
time_labels = [txt for _ = 1:Nrows, txt ∈ [L"B_\mathrm{sim}", L"B_{\beta}",
    L"B_{\mathcal{F}}", L"B_\mathrm{ff}",
    L"B_\mathrm{dyn ↓}", L"B_\mathrm{dyn ↑}"]
]

log_map = falses(Nrows * Ncols)
annotate_scale = trues(Nrows * Ncols)

scale_kpc = 1_000.0
scale_label = "1 Mpc"

plot_name = plot_path * "Fig5a.pdf"

# plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
#     vmin_arr, vmax_arr, plot_name,
#     cutoffs=[vmin_arr[1] for i = 1:Ncols],
#     #mask_bad=trues(Nrows * Ncols),
#     upscale=2.5,#0.6,
#     cb_label_offset=0.6,
#     dpi=400,
#     transparent=false,
#     ticks_color="k",
#     annotation_color="k",
#     colorbar_mode="single",
#     read_mode=4;
#     time_labels, annotate_time,
#     log_map,
#     annotate_scale,
#     scale_kpc,
#     scale_label,
# )


"""
    Radial Profiles
"""
function plot_radial_B(plot_name)

    Bfield_models = [L"B_\mathrm{sim}",
        L"B_{\beta}",
        L"B_{\mathcal{F}}",
        L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}",
        L"B_\mathrm{dyn ↑}"]

    Bfield_names = ["B_sim", "B_beta", "B_vturb", "B_FF", "B_dyn_l", "B_dyn_h"]
    folders = "radial_B/" .* ["box", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_inj"]

    sim_names = ["SLOW-CR3072", "Coma", "Coma-" * L"D_\mathrm{pp}" * "-Low",
        "Coma-" * L"D_\mathrm{pp}" * "-High"]

    bonafede = readdlm(data_path * "bonafede_sorted.dat")

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

        if col <= 2

            for i_sim = 1:length(sim_names)

                fi = data_path * "/" * folders[i_sim] * "/" * Bfield_names[col] * ".dat"
                d = readdlm(fi)
                plot(d[:,1], 1.e6 .* d[:,2], 
                    alpha=alpha_ref,
                    color=sm2.to_rgba(i_sim), lw=lw,
                    label=sim_names[i_sim])

            end
        end


        if 3 <= col <= 4

            for i_sim = 2:2#length(sim_names)

                fi = data_path * "/" * folders[i_sim] * "/" * Bfield_names[col] * ".dat"
                d = readdlm(fi)
                plot(d[:,1], 1.e6 .* d[:,2], 
                    alpha=alpha_ref,
                    color=sm2.to_rgba(i_sim), lw=lw,
                    linestyle="--")

            end
        end

        for i_sim = 1:length(sim_names)

            fi = data_path * "/" * folders[i_sim] * "/" * Bfield_names[col] * "_rescaled.dat"

            if !isfile(fi)
                continue
            end

            d = readdlm(fi)
            plot(d[:,1], 1.e6 .* d[:,2], 
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

folders = ["box", "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_HB07"]
plot_name = plot_path * "Fig05b.pdf"

plot_radial_B(plot_name)

