using GadgetIO, GadgetUnits
using Printf
using ProgressMeter
using DelimitedFiles
using Base.Threads
using Statistics
using Distributions
using StatsBase


const global h = read_header("/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20e/snapdir_074/snap_074.0")
const global GU = GadgetPhysical(h)


include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))


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

    data_path = "/e/ocean2/users/lboess/PaperRepos/2310.13734/data/B_histograms/"

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

sim_names = ["box", "zoom_inj", "zoom_dpp"]

snap_bases = ["/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000",
            "/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20e/snapdir_074/snap_074",
            "/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20eDpp/snapdir_074/snap_074"]

use_keys = [true, false, false]

Bfield_functions = [#Bfield_sim, Bfield_FF, Bfield_Beta, 
                    Bfield_vturb#, Bfield_dyn_l, Bfield_dyn_h
                    ]
Bfield_names = [#"sim", "ff", "beta", 
                "vturb1"#, "dyn_l", "dyn_h"
                ]

for i in 1:length(snap_bases)
    for j = 1:length(Bfield_functions)
        write_histograms(snap_bases[i], sim_names[i] * "_" * Bfield_names[j], use_keys[i], Bfield_functions[j])
    end
end