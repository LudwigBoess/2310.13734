include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using GadgetIO
using PyPlot, PyPlotUtility
using PyCall

cmr = pyimport("cmasher")
PE = pyimport("matplotlib.patheffects")
@info "done"

"""
    plot_text_clusters(ax, clusters)

Plot text annotations for clusters in the given list of `clusters` on the current plot axis.
"""
function plot_text_clusters(ax, clusters)

    for cl ∈ clusters
        txt = ax.text(cl.i, cl.j, cl.name, fontsize=20, color="w")
        txt.set_path_effects([PE.withStroke(linewidth=2, foreground="k")])
    end

    return ax
end

"""
    ClusterAnnotation

Struct to store the name and pixel coordinates of a cluster.
"""
struct ClusterAnnotation
    name::String
    xpix::Int64
    ypix::Int64
end


cluster = [ClusterAnnotation("A1644", 2560 - 250, 1479 - 30),
    ClusterAnnotation("A85", 3787 - 2448, 194 - 80),
    ClusterAnnotation("A119", 3678 - 2048 - 950, 266),
    ClusterAnnotation("A347", 285 + 30, 967 - 30),
    ClusterAnnotation("A496", 3039 - 50, 306 - 100),
    ClusterAnnotation("A539", 3937 + 20, 1024),
    ClusterAnnotation("A576", 315 + 10, 1267 + 50),
    ClusterAnnotation("A644", 3454 + 30, 1201 - 20),
    ClusterAnnotation("A754", 3454 - 70, 1201 + 100),
    ClusterAnnotation("A1185", 3089, 1876 - 80),
    ClusterAnnotation("A1367", 2803 + 20, 1896 + 20),
    ClusterAnnotation("A3158", 3026 + 20, 407 - 30),
    ClusterAnnotation("A1795", 2460 - 650, 2039 - 80),
    ClusterAnnotation("A2029", 1144 + 850, 1634),
    ClusterAnnotation("A2063", 1274 + 600, 1541 - 50),
    ClusterAnnotation("A2065", 283 + 1300, 1789 + 30),
    ClusterAnnotation("A2107", 1007 + 570, 1631),
    ClusterAnnotation("A2199", 1231 + 300, 1474 - 90),
    ClusterAnnotation("A2256", 441 + 300, 1540 - 30),
    ClusterAnnotation("A1736", 2639 - 250, 1359 - 50),
    ClusterAnnotation("A2593", 22 + 1300, 289 - 30),
    ClusterAnnotation("A2634", 553 + 50, 629),
    ClusterAnnotation("A2665", 3918 - 3000, 526 - 30),
    ClusterAnnotation("A2734", 3636 - 1900, 92 + 50),
    ClusterAnnotation("A2877", 3233 - 1130, 45),
    ClusterAnnotation("A3266", 3074 - 250, 493 - 30),
    ClusterAnnotation("A3571", 2688 - 80, 1285 - 80),
    ClusterAnnotation("A3376", 3199 + 20, 535 - 20),
    ClusterAnnotation("A3391", 3151 + 20, 662 + 20),
    ClusterAnnotation("A3532", 2666 + 20, 1430 - 30),
    ClusterAnnotation("A3581", 2684 + 20, 1484 + 20),
    ClusterAnnotation("A3667", 2181 + 20, 482),
    ClusterAnnotation("AWM7", 238 - 50, 910 - 90),
    ClusterAnnotation("Centaurus", 2869 + 50, 1299 + 50),
    ClusterAnnotation("Coma", 2586, 1970 + 25),
    ClusterAnnotation("Fornax", 3615 + 30, 707 + 30),
    ClusterAnnotation("Hydra", 3013 + 50, 1229 + 30),
    ClusterAnnotation("Norma", 2613 + 30, 1007),
    ClusterAnnotation("Perseus", 229 - 80, 1017 + 50),
    ClusterAnnotation("Virgo", 2567 - 200, 1867 - 130),
    ClusterAnnotation("A2572a", 61 + 1050, 358),
    ClusterAnnotation("2A0335+096", 295, 444)]


dpi = 100
file_ending = "pdf"
contour_file = map_path * "allsky/cluster_contours_slow_1_gal.fits"


filename = "CReP_seed"
im_cmap = cmr.freeze
cb_label = "CR Electron Energy Density  " * L"\epsilon_{\mathrm{CR,e} \in 5-500 \mathrm{ MeV}}" * " [erg cm" * L"^{-3}" * "]"
clim_arr = [1.e-22, 1.e-17]
plot_name = "fig01a"

filename = map_path * "allsky/allsky_" * filename * "_slow_1_gal.fits"
plot_name = plot_path * "$plot_name.$file_ending"

plot_single_allsky(filename, im_cmap, cb_label, clim_arr, plot_name,
    time_label="",
    log_map=true,
    Npixels=2048,
    contour_linestyle="-",
    contour_alpha=1.0,
    contour_color="w",
    annotations=cluster,
    dpi=dpi,
    cutoff=clim_arr[1],
    origin="lower"
    ;
    contour_file
)


filename = "CReE_40Mpc"
im_cmap = cmr.cosmic
cb_label = "CR Electron Energy " * L"E_{\mathrm{CR},e > 1 GeV}" * " [erg cm" * L"^{-2}" * "]"
clim_arr = [1.e2, 1.e6]
plot_name = "Fig01b"

filename = map_path * "allsky/allsky_" * filename * "_slow_1_gal.fits"
plot_name = plot_path * "$plot_name.$file_ending"

plot_single_allsky(filename, im_cmap, cb_label, clim_arr, plot_name,
    time_label="",
    log_map=true,
    Npixels=2048,
    contour_linestyle="-",
    contour_alpha=1.0,
    contour_color="w",
    dpi=dpi,
    cutoff=clim_arr[1],
    origin="lower"
)