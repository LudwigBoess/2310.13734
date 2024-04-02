using GadgetIO, GadgetUnits
using SpectralCRsUtility
using Printf
using ProgressMeter
using DelimitedFiles
using Base.Threads
using Statistics
using Distributions
using StatsBase
using PyPlot, PyPlotUtility


const global h = read_header("/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/mhd_cr8p24eDpp_5e-17/snapdir_012/snap_012.0")
const global GU = GadgetPhysical(h)


include(joinpath(@__DIR__, "allsky", "Bfld.jl"))


function read_data(snap_base, use_keys)

    blocks = ["POS", "MASS", "HSML", "RHO", "U", 
              "BFLD", "VRMS", "CReP"]

    # found by selecting particles by hand
    cylinder = GadgetCylinder([247.75, 332.0, 248.5] .* 1.e3, 
                              [248.5, 332.75, 249.25] .* 1.e3, 
                              2.5e3)

    d = read_particles_in_geometry(snap_base, blocks, cylinder; use_keys)

    sel = findall(d["CReP"] .> 0.0)

    # iterate blocks of dict to remove filtered out particles
    for block ∈ keys(d)
        # get all dimensions but the last (the last is masked)
        colons = repeat([:], ndims(d[block]) - 1)

        d[block] = d[block][colons..., sel]
    end

    return d
end

function get_B(data, Bfield_function)

    B = Vector{Float64}(undef, size(data["POS"], 2))

    @threads for i = 1:size(data["POS"], 2)
        B[i] = log10(Bfield_function(data, i))
    end

    return B
end

"""
    get_histograms(bins, Q)

Computes the histograms for `M_s` of all shocked particles
"""
function get_histograms(bins, Q)

    # fit histograms
    hist = fit(Histogram, Q, bins)

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

"""
    write_histograms(snap, Nbins=50)

Constructs the B number histograms of all shocked particles and writes them to a txt file.
"""
function write_histograms(snap_base, sim_name, use_keys, Bfield_function, Nbins=50)

    data_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/B_histograms/"

    # construct boundaries and centers for binning
    B_bins, B_centers = construct_bins_and_centers(-12, -4, Nbins)

    # reads the simulation data
    data = read_data(snap_base, use_keys)

    B = get_B(data, Bfield_function)

    # construct histograms
    B_binned = get_histograms(B_bins, B)

    # write into single array for saving
    d = [ B_centers B_binned ]

    # save as txt file
    data_file = data_path * "$(sim_name)_B_histograms.txt"
    writedlm(data_file, d)
end

snap_bases = ["/gpfs/work/pn68va/di67meg/LocalUniverse/snapdir_036/snap_036",
    "/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/mhd_cr8p24e/snapdir_012/snap_012",
    "/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/mhd_cr8p24eDpp_1e-17/snapdir_074/snap_074",
    "/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/mhd_cr8p24eDpp_5e-17/snapdir_012/snap_012"]

sim_names = ["box", "zoom_inj", "zoom_dpp_low", "zoom_dpp_high"]
use_keys = [true, false, false, false]
Bfield_functions = [Bfield_sim, Bfield_FF, Bfield_Beta, Bfield_vturb, Bfield_dyn_l, Bfield_dyn_h]
Bfield_names = ["sim", "ff", "beta", "vturb", "dyn_l", "dyn_h"]

# for i in 1:length(snap_bases)
#     for j = 5:5#1:length(Bfield_functions)
#         write_histograms(snap_bases[i], sim_names[i] * "_" * Bfield_names[j], use_keys[i], Bfield_functions[j])
#     end
# end


function read_histograms(data_file)
    d = readdlm(data_file)
    return d[:, 1], d[:, 2]./sum(d[:, 2])
end

function plot_histograms(sim_names, Bfield_names, plot_name)

    data_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/B_histograms/"

    Bfield_models = [L"B_\mathrm{sim}",
                    L"B_{\beta = 50}",
                    L"B_{\mathcal{F} = 0.1}",
                    L"B_\mathrm{ff}",
                    L"B_\mathrm{dyn ↓}",
                    L"B_\mathrm{dyn ↑}"]

    sim_labels = ["SLOW-CR3072" * L"^3",
                    "Coma",
                    "Coma-" * L"D_{\mathrm{pp}}" * "-Low",
                    "Coma-" * L"D_{\mathrm{pp}}" * "-High"]

    Nrows = length(sim_names)

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=6.5))
    sm2.set_array([])

    x_pixels = 1500
    fig = get_figure(0.25; x_pixels)
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

        axvline(log10(3.2e-6), color="k", linestyle="--", label=L"B_\mathrm{CMB}")
        text(-11.5, 1.0, sim_labels[i+1], fontsize=16, horizontalalignment="left")

        for j = 1:length(Bfield_models)
            sim_name = sim_names[i+1] * "_" * Bfield_names[j]
            data_file = data_path * "$(sim_name)_B_histograms.txt"
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


sim_names = ["box", "zoom_inj", "zoom_dpp_low", "zoom_dpp_high"]
Bfield_names = ["sim", "ff", "beta", "vturb", "dyn_l", "dyn_h"]
plot_name = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/B_histograms.pdf"

plot_histograms(sim_names, Bfield_names, plot_name)