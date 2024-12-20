include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "shared.jl"))

using SPHtoGrid
using Printf
using GadgetIO, GadgetUnits
using Unitful, UnitfulAstro
using PyPlot, PyPlotUtility
using PyCall
cm = pyimport("cmasher")
using ColorSchemes


"""
    write_baseline_image(map_path, snap, Bname, smooth_size)


"""
function write_baseline_image(map_path, snap, Bname, smooth_size)

    filename = "synch_F_beam_1'_144MHz" * "_$Bname"

    file = map_path * "coma_20Mpc_$snap.$filename.xz.fits"

    println("loading data")
    im_default, par, snap_num, units = read_fits_image(file)
    im_smoothed = PyPlotUtility.smooth_map!(im_default, [smooth_size, smooth_size], par)

    println("constructing baseline image")
    im_baseline = (im_default .- im_smoothed) ./ im_default

    filename_out = map_path * "coma_20Mpc_$snap.baseline_$Bname.xz.fits"
    println("saving baseline image: $filename_out")
    write_fits_image(filename_out, im_baseline,
        par,
        units="[]",
        snap=snap_num)

    println("constructing baseline image")
    im_baseline = (im_default .- im_smoothed) ./ im_default

    filename_out = map_path * "coma_20Mpc_$snap.synch_F_beam_1'_144MHz_boosted_$Bname.xz.fits"

    im_boosted = @. im_baseline * im_default + im_default
    println("saving baseline image: $filename_out")
    write_fits_image(filename_out, im_boosted,
        par,
        units="[]",
        snap=snap_num)
    println("done")

end

# Bnames = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]
# folders = ["box", 
#         "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_HB07"] .* "/"
# snaps = ["036", 
#     "012", "074", "012", "012"]

# snapshot_file = "/path/to/snapshot/snap_000"

# # read the Gadget file header
# h = read_header(snapshot_file)

# # define a cosmology from the header
# c = cosmology(h)

# # observed redshift of Coma
# z = 0.0231

# # lofar
# nu = 144.0u"MHz"
# Dmin = 0.1u"km"
# Smax = 0.6u"c0"/nu/Dmin |> upreferred
# theta = Smax * 1.0u"rad" |> u"arcminute"

# # scale of missing baseline contribution in kpc
# smooth_scale = arcmin_to_kpc(c, theta, z)

# println("Lofar missing baseline scale in Coma: ", theta, "  ", smooth_scale)


# uncomment to run baseline computation
# for i ∈ 4:4, Bname ∈ Bnames 
#     println(folders[i], " ", Bname)
#     write_baseline_image(map_path * folders[i], snaps[i], Bname, smooth_scale)
# end


Bfield_names = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]
files = @. map_path * "coma/zoom_inj/coma_20Mpc_$(@sprintf("%03i", 12)).baseline_" * Bfield_names * ".xz.fits"

Ncols = length(Bfield_names)
Nrows = 1

vmin_arr = [-1.0]
vmax_arr = [1.0]

im_cmap = [cm.prinsenvlag]

cb_labels = ["Relative Missing Baseline Contribution  " * L"(I-I')/I"]

annotate_time = trues(Nrows * Ncols)
time_labels = [txt for _ = 1:Nrows, txt ∈ [L"B_\mathrm{sim}", L"B_{\beta}",
    L"B_{\mathcal{F}}", L"B_\mathrm{ff}",
    L"B_\mathrm{dyn ↓}", L"B_\mathrm{dyn ↑}"]
]

log_map = falses(Nrows * Ncols)
annotate_scale = trues(Nrows * Ncols)

scale_kpc = 5_000.0
scale_label = "5 Mpc"

plot_name = plot_path * "Fig15a.pdf"

plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
    vmin_arr, vmax_arr, plot_name,
    cutoffs=[vmin_arr[1] for i = 1:Ncols],
    #mask_bad=trues(Nrows * Ncols),
    upscale=3.0,
    cb_label_offset=0.6,
    dpi=400,
    transparent=false,
    ticks_color="k",
    colorbar_mode="single",
    annotation_color="k",
    read_mode=1,
    mask_bad=trues(Nrows * Ncols);
    time_labels, annotate_time,
    log_map,
    annotate_scale,
    scale_kpc,
    scale_label
)


"""
    Boosted
"""
Bfield_names = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]
files = @. map_path * "coma/zoom_inj/coma_20Mpc_$(@sprintf("%03i", 12)).synch_F_beam_1'_144MHz_boosted_" * Bfield_names * ".xz.fits"

Ncols = length(Bfield_names)
Nrows = 1

    vmin_arr = [1.e-8]
    vmax_arr = [1.e2]

    bk = ColorMap("bukavo", ColorSchemes.bukavu[1:end])

    im_cmap = [bk]

    cb_labels = ["Synchrotron Flux + Response " * L"F_{ν,144 \mathrm{ MHz}}" *
                 " [mJy arcmin" * L"^{-2}" * "]"]

annotate_time = trues(Nrows * Ncols)
time_labels = [txt for _ = 1:Nrows, txt ∈ [L"B_\mathrm{sim}", L"B_{\beta}",
    L"B_{\mathcal{F}}", L"B_\mathrm{ff}",
    L"B_\mathrm{dyn ↓}", L"B_\mathrm{dyn ↑}"]
]

log_map = trues(Nrows * Ncols)
annotate_scale = trues(Nrows * Ncols)

scale_kpc = 5_000.0
scale_label = "5 Mpc"

plot_name = plot_path * "Fig15b.pdf"

plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
    vmin_arr, vmax_arr, plot_name,
    cutoffs=[vmin_arr[1] for i = 1:Ncols],
    #mask_bad=trues(Nrows * Ncols),
    upscale=3.0,
    cb_label_offset=0.6,
    dpi=400,
    transparent=false,
    ticks_color="k",
    colorbar_mode="single",
    read_mode=1,
    mask_bad=falses(Nrows * Ncols);
    time_labels, annotate_time,
    log_map,
    annotate_scale,
    scale_kpc,
    scale_label
)