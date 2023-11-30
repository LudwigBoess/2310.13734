include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using PyCall
using ColorSchemes
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
#folders = ["box", "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "analytic_HB"]
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
#folders = ["coma_filament", "zoom_dpp_1e-17", "zoom_dpp_1e-16"]
plot_10Mpc_slope(map_path *  "coma/", folders)