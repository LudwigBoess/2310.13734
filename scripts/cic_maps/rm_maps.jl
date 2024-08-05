using GadgetIO, GadgetUnits
using Printf 
using PyPlot, PyPlotUtility
using PyCall
cm = pyimport("cmasher")


plot_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/"

# Bfield_names = ["B_sim", "B_beta", "B_vturb",
#                 "B_FF", "B_dyn_l", "B_dyn_h"]

Bfield_names = ["B_sim", "B_beta", 
    "B_vturb_rescaled",
    "B_FF_rescaled", "B_dyn_l_rescaled", "B_dyn_h_rescaled"]


files = @. "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/coma/zoom_dpp_1e-17/coma_5Mpc_$(@sprintf("%03i", 74)).RM_" * Bfield_names * ".a.y.fits"

Ncols = length(Bfield_names)
Nrows = 1

vmin_arr = [-300.0]
vmax_arr = [300.0]

im_cmap = ["seismic"] #cm.copper]

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

plot_name = plot_path * "FigRM.pdf"

plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
    vmin_arr, vmax_arr, plot_name,
    cutoffs=[vmin_arr[1] for i = 1:Ncols],
    #mask_bad=trues(Nrows * Ncols),
    upscale=2.5,#0.6,
    cb_label_offset=0.6,
    dpi=400,
    transparent=false,
    ticks_color="k",
    annotation_color="k",
    colorbar_mode="single",
    read_mode=4;
    time_labels, annotate_time,
    log_map,
    annotate_scale,
    scale_kpc,
    scale_label,
)
