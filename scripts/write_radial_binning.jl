using GadgetIO, GadgetUnits
using Printf
using DelimitedFiles
using ProgressMeter

include(joinpath(@__DIR__, "allsky", "Bfld.jl"))
include(joinpath(@__DIR__, "phase_maps", "bin_1D.jl"))

function get_r(data, GU)
    
    pot_min = findmin(data["POT"])[2]
    center = data["POS"][:,pot_min]
    r = Vector{Float64}(undef, length(data["POS"][1,:]))

    @threads for i = 1:length(r)
        r[i] = 0.0
        for j = 1:3
            r[i] += (data["POS"][j,i] - center[j])^2
        end  
        r[i] = sqrt(r[i])
    end

    return r .* GU.x_physical
end


function get_bin_centers(r_lim, Nbins)

    # get logarithmic bin spacing
    dbin = (log10(r_lim[2]) - log10(r_lim[1])) / Nbins

    # define bin centers 
    return [r_lim[1] * 10^((i + 1 / 2) * dbin) for i = 0:Nbins-1]
end


function get_B_binned(r, B)
    Nbins = 50
    r_lim = [30, 1.e4]

    count, B_bin = bin_1D_log(r, r_lim, B, calc_mean=true, calc_sigma=false, show_progress=true; Nbins)

    return get_bin_centers(r_lim, Nbins), B_bin
end


function write_binning(filename, data, Bfield_function)

    r = get_r(data, GU)
    B = Vector{Float64}(undef, length(r))

    @threads for i = 1:length(r)
        B[i] = Bfield_function(data, i)
    end

    r_centers, B_bin = get_B_binned(r, B)

    writedlm(filename, [r_centers B_bin])
end

function read_data(snap_base, use_keys)

    blocks = ["POS", "VEL", "HSML", "RHO", "U", "MASS", "VRMS",
        "BFLD", "POT"]

    gpos = [245040.45, 327781.84, 246168.69]
    rvir = 2177.5625

    read_particles_in_volume(snap_base, blocks, gpos, 2rvir, 
                            use_keys=use_keys)
end

folders = ["box", "zoom_inj", "zoom_dpp"]
data_paths = "/e/ocean2/users/lboess/PaperRepos/2310.13734/data/radial_B/" .* folders .* "/"

snap_base = ["/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000",
            "/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20e/snapdir_074/snap_074",
            "/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20eDpp/snapdir_074/snap_074"]

use_keys = [true, false, false]

const global GU = GadgetPhysical(read_header("/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20e/snapdir_074/snap_074.0"))

for j = 1:length(data_paths)

    data = read_data(snap_base[j], use_keys[j])

    Bfield_functions = [#Bfield_sim, Bfield_Beta, 
                        Bfield_vturb#, Bfield_FF, Bfield_dyn_l, Bfield_dyn_h
                        ]
    Bfield_names = [#"B_sim", "B_beta", 
                    "B_vturb1"#, "B_FF", "B_dyn_l", "B_dyn_h"
                    ]

    for i = 1:length(Bfield_functions)
        write_binning(data_paths[j] * Bfield_names[i] * ".dat", data, Bfield_functions[i])
    end

end