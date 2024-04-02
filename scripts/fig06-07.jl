include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using PyCall
using ColorSchemes
using Distributions
using StatsBase
cm = pyimport("cmasher")
@info "done"

function plot_10Mpc_col(map_path, folders)

    filenames = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]

    snaps = ["036", "012", "074", "012", "074"]

    Ncols = length(filenames)
    Nrows = 5

    files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_F_beam_1'_$filename.xz.fits"
             for i ∈ 1:Nrows, filename ∈ filenames]


    vmin_arr = [1.e-8]
    vmax_arr = [1.e2]

    bk = ColorMap("bukavo", ColorSchemes.bukavu[1:end])

    im_cmap = [bk]

    cb_labels = ["Synchrotron Flux  " * L"F_{ν,144 \mathrm{ MHz}}" *
                 " [mJy arcmin" * L"^{-2}" * "]"]


    annotate_time = trues(Nrows * Ncols)
    time_labels = [txt for _ = 1:Nrows, txt ∈ [L"B_\mathrm{sim}", L"B_{\beta = 50}",
        L"B_{\mathcal{F} = 0.1}", L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}", L"B_\mathrm{dyn ↑}"]
    ]

    time_labels[1] = L"B_\mathrm{sim}" * "\n" * "SLOW-CR3072" * L"^3"
    time_labels[2] = "Coma"
    time_labels[3] = "Coma-" * L"D_\mathrm{pp}" * "-Low"
    time_labels[4] = "Coma-" * L"D_\mathrm{pp}" * "-High"
    time_labels[5] = "H&B (2007)"

    log_map = trues(Nrows * Ncols)
    annotate_scale = trues(Nrows * Ncols)

    scale_kpc = 5_000.0
    scale_label = "5 Mpc"

    plot_name = plot_path * "Fig05.pdf"

    plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
        vmin_arr, vmax_arr, plot_name,
        cutoffs=[vmin_arr[1] for i = 1:Ncols],
        mask_bad=trues(Nrows * Ncols),
        upscale=0.6,
        cb_label_offset=0.6,
        dpi=400,
        transparent=false,
        ticks_color="k",
        colorbar_mode="single";
        time_labels, annotate_time,
        log_map,
        annotate_scale,
        scale_kpc,
        scale_label,
    )
    GC.gc()
end


map_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/"
folders = ["box", "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_HB07"]
plot_10Mpc_col(map_path * "coma/", folders)


function plot_10Mpc_slope(map_path, folders)

    filenames = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]

    snaps = ["036", "012", "074", "012", "074"]

    Ncols = length(filenames)
    Nrows = 5

    files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_slope_$filename.xz.fits"
             for i ∈ 1:Nrows, filename ∈ filenames]

    vmin_arr = [-3.0]
    vmax_arr = [0.0]

    im_cmap = ["Spectral_r"]

    cb_labels = ["Synchrotron Spectral Slope  " * L"\alpha_\mathrm{144 MHz}^\mathrm{944 MHz}"]


    annotate_time = trues(Nrows * Ncols)
    time_labels = [txt for _ = 1:Nrows, txt ∈ [L"B_\mathrm{sim}", L"B_{\beta = 50}",
        L"B_{\mathcal{F} = 0.1}", L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}", L"B_\mathrm{dyn ↑}"]
    ]

    time_labels[1] = L"B_\mathrm{sim}" * "\n" * "SLOW-CR3072" * L"^3"
    time_labels[2] = "Coma"
    time_labels[3] = "Coma-" * L"D_\mathrm{pp}" * "-Low"
    time_labels[4] = "Coma-" * L"D_\mathrm{pp}" * "-High"
    time_labels[5] = "H&B (2007)"

    log_map = falses(Nrows * Ncols)
    annotate_scale = trues(Nrows * Ncols)

    scale_kpc = 5_000.0
    scale_label = "5 Mpc"
    circle_color = "gray"

    plot_name = plot_path * "Fig06.pdf"

    plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
        vmin_arr, vmax_arr, plot_name,
        cutoffs=[vmin_arr[1] for i = 1:Ncols],
        mask_bad=trues(Nrows * Ncols),
        upscale=0.6,
        cb_label_offset=0.6,
        dpi=400,
        transparent=false,
        ticks_color="k",
        colorbar_mode="single",
        annotation_color="k";
        time_labels, annotate_time,
        log_map,
        annotate_scale,
        scale_kpc,
        scale_label,
        circle_color
    )
    GC.gc()
end


map_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/"
folders = ["box", "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_HB07"]
plot_10Mpc_slope(map_path *  "coma/", folders)


"""
    get_histograms(bins, Q)

Computes the histograms for `M_s` of all shocked particles
"""
function get_histograms(bins, Q, weights)

    # fit histograms
    hist = fit(Histogram, Q, weights, bins)

    return hist.weights
end

"""
    construct_bins_and_centers(bin_min, bin_max, Nbins)

Returns the bins and their center points.
"""
function construct_bins_and_centers(bin_min, bin_max, Nbins)

    # construct bin boundaries
    bins = LinRange(bin_min, bin_max, Nbins + 1)

    # construct bin centers
    bin_centers = Vector{Float64}(undef, Nbins)

    for i = 1:Nbins
        bin_centers[i] = 0.5 * (bins[i] + bins[i+1])
    end

    bins, bin_centers
end

function read_data(α_file, jν_file)

    bins, bin_centers = construct_bins_and_centers(-3.0, 0.0, 40)

    map, par, snap_num, units = PyPlotUtility.read_fits_image(α_file)
    w_map, par, snap_num, units = PyPlotUtility.read_fits_image(jν_file)

    sel = findall(.!isnan.(map))

    println("pixels above threshold: ", length(findall(map[sel] .> -0.5)))

    if length(findall(map[sel] .> -0.5)) > 0
        println(map[findall(map[sel] .> -0.5)])
    end

    w = AnalyticWeights(ones(length(sel)))
    hist = get_histograms(bins, map[sel], w)


    sum_jnu = sum(w_map[sel])
    w = AnalyticWeights(w_map[sel] ./ sum_jnu)
    hist_w = get_histograms(bins, map[sel], w)

    return bin_centers, hist, hist_w
end

function plot_slope_histograms(map_path, folders, plot_name)

    Bfield_models = [L"B_\mathrm{sim}",
        L"B_{\beta = 50}",
        L"B_{\mathcal{F} = 0.1}",
        L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}",
        L"B_\mathrm{dyn ↑}"]

    snaps = ["036", "012", "074", "012", "074"]


    Bfield_filenames = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]


    sim_names = ["SLOW-CR3072", "Coma", "Coma-" * L"D_\mathrm{pp}" * "-Low",
        "Coma-" * L"D_\mathrm{pp}" * "-High", "H&B (2007)"]

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=4.5))
    sm2.set_array([])

    alpha_ref = 1.0
    lw = 5.0
    x_pixels = 600
    legend_font_size = 30
    fig = get_figure(5.5; x_pixels)
    plot_styling!(x_pixels; legend_font_size)

    gs = plt.GridSpec(1, 6, figure=fig)

    for col = 1:length(Bfield_models)

        println("col $col")
        subplot(get_gs(gs, 0, col-1))
        ax = gca()

        axis_ticks_styling!(ax)
        ax.set_yscale("log")
        ax.set_xlim([-3.0, 0.0])
        ax.set_ylim([1.e-3, 1.0])
        ax.set_facecolor("white")

        ax.set_xticklabels([-3.0, -2.5, -2.0, -1.5, -1.0, -0.5])

        if col == 3
            xlabel("Synchrotron Spectral Slope  " * L"\alpha_\mathrm{144 MHz}^\mathrm{944 MHz}",
                    fontsize=legend_font_size)
            ax.xaxis.set_label_coords(0.8, -0.05)
        end
        if col == 1
            ylabel("Rel. Number of Pixels  " * L"N/N_\mathrm{tot}", fontsize=legend_font_size)
            #ax.yaxis.set_label_coords(-0.2, 0.05)
        else
            ax.set_yticklabels([])
        end

        for i_sim = 1:length(sim_names)-1
            α_file = map_path * "$(folders[i_sim])/coma_20Mpc_$(snaps[i_sim]).synch_slope_$(Bfield_filenames[col]).xz.fits"
            jν_file = map_path * "$(folders[i_sim])/coma_20Mpc_$(snaps[i_sim]).synch_F_beam_1'_$(Bfield_filenames[col]).xz.fits"
            
            bins, hist, hist_w = read_data(α_file, jν_file)
            hist ./= sum(hist)
            plot(bins, hist, alpha=alpha_ref,
                color=sm2.to_rgba(i_sim), lw=lw, 
                label=sim_names[i_sim])

            plot(bins, hist_w, alpha=alpha_ref,
                color=sm2.to_rgba(i_sim), lw=lw,
                linestyle=":")
        end

        α_file  = map_path * "$(folders[end])/coma_20Mpc_$(snaps[end]).synch_slope_$(Bfield_filenames[col]).xz.fits"
        jν_file = map_path * "$(folders[end])/coma_20Mpc_$(snaps[end]).synch_F_beam_1'_$(Bfield_filenames[col]).xz.fits"

        bins, hist, hist_w = read_data(α_file, jν_file)
        hist ./= sum(hist)

        plot(bins, hist, alpha=alpha_ref,
            color="gray", lw=lw, 
            label=sim_names[end])

        plot(bins, hist_w, alpha=alpha_ref,
            color="gray", lw=lw,
            linestyle=":")

        ax.text(-2.8, 0.5, Bfield_models[col], fontsize=legend_font_size,
            horizontalalignment="left", verticalalignment="center")

        ax.axvspan(-0.5, 0.0, color="k", alpha=0.2)

    end

    plot([0.0], [0.0], color="k", lw=lw, label=L"\sum \: N")
    plot([0.0], [0.0], color="k", lw=lw, linestyle=":", label=L"\sum \: j_\nu")

    handles, labels = gca().get_legend_handles_labels()
    order = [1, 6, 2, 7, 3, 4, 5]
    l = legend([handles[idx] for idx in order], [labels[idx] for idx in order],
        frameon=false, bbox_to_anchor=(-2.0, -0.5), ncol=5, loc="lower center")

    subplots_adjust(hspace=0.0, wspace=0.0)
    savefig(plot_name, bbox_inches="tight", transparent=false)
    close(fig)
end

map_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/coma/"
folders = ["box", "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_HB07"]
plot_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/"
plot_name = plot_path * "Fig06b.pdf"

plot_slope_histograms(map_path, folders, plot_name)

filenames = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]

snaps = ["036", "012", "074", "012", "074"]

Ncols = length(filenames)
Nrows = 5

α_files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_slope_$filename.xz.fits"
         for i ∈ 1:Nrows, filename ∈ filenames]

jν_files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_F_beam_1'_$filename.xz.fits"
             for i ∈ 1:Nrows, filename ∈ filenames]


