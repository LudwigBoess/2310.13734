include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using GadgetIO
using PyPlot, PyPlotUtility
using PyCall
using Unitful, UnitfulAstro
@info "done"

# plot settings
dpi = 100
file_ending = "png"

filenames = map_path * "allsky/" .* ["allsky_B_sim_slow_1_gal.fits",
    "allsky_B_beta50_slow_1_gal.fits",
    "allsky_B_01Pturb_slow_1_gal.fits",
    "allsky_B_FF_slow_1_gal.fits",
    "allsky_B_dyn_l_slow_1_gal.fits",
    "allsky_B_dyn_h_slow_1_gal.fits"]

time_label = [L"B_\mathrm{sim}", L"B_{\beta}",
    L"B_{\mathcal{F}}", L"B_\mathrm{ff}",
    L"B_\mathrm{dyn ↓}", L"B_\mathrm{dyn ↑}"]

im_cmap = "magma"
cb_label = "Magnetic Field Strength  " * L"\vert B \vert" * " [G]"
clim_arr = [1.e-9, 1.0e-5]
plot_name = plot_path * "Fig02.$file_ending"
plot_multiple_allsky(filenames, im_cmap, cb_label, clim_arr, plot_name,
                    time_labels=time_label,
                    log_map=true,
                    Npixels=2048,
                    dpi=dpi,
                    origin="lower",
                    ticks_color="w",
                    Ncols=3,
                    Nrows=2,
                    transparent=true
                )
