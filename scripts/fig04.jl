include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using PyCall
using ColorSchemes
cm = pyimport("cmasher")
@info "done"

function plot_general_quantities(map_path, folders)

    filenames = ["rho", "B", "CReE_gt1GeV_L", "D0_max"]

    snaps = ["036", "012"]

    Ncols = length(filenames)
    Nrows = 2

    files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).$filename.xz.fits"
             for i ∈ 1:2, filename ∈ filenames]


    im_cmap = ["viridis", "magma", 
        cm.cosmic, "gnuplot2"]

    cb_labels = [L"\Sigma_g" * " [cm" * L"^{-2}" * "]",
        L"B" * " [" * L"\mu" * "G]",
        L"E_{\mathrm{CR},e > 1 GeV}" * " [erg]",
        L"D_0" * " [s" * L"^{-1}" * "]"]


    annotate_time = trues(Nrows * Ncols)
    time_labels = ["" for _ = 1:Nrows*Ncols]
    time_labels[1] = "SLOW-CR3072" * L"^3"
    time_labels[2] = "Coma"

    log_map = trues(Ncols * Nrows)
    annotate_scale = trues(Ncols * Nrows)
    #r_circles = [2642.8608, 1837.9205]
    r_circles = [3895.850659851734, 2709.2852602308144]

    vmin_arr = [1.e-4, 2.e-4, 8.e48, 4e-18]
    vmax_arr = [1.e-1, 8.0, 4.e53, 1.e-15]

    plot_name = plot_path * "fig04.pdf"

    plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
        vmin_arr, vmax_arr, plot_name, upscale=1.2,
        cutoffs=vmin_arr,
        mask_bad=trues(Ncols),
        scale_kpc=5_000.0,
        scale_label="5 Mpc",
        cb_label_offset=1.1,
        dpi=400,
        circle_alpha=0.8;
        time_labels, annotate_time,
        log_map,
        r_circles,
        annotate_scale
    )

end

folders = ["box", "zoom_inj"]
plot_general_quantities(map_path * "coma/", folders)