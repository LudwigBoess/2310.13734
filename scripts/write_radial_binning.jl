using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using DelimitedFiles

include(joinpath(@__DIR__, "allsky", "Bfld.jl"))

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

function read_data(sim_path, snap, use_keys)

    snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
    blocks = ["POS", "VEL", "HSML", "RHO", "U", "MASS", "VRMS",
        "BFLD", "POT"]

    gpos = [245040.45, 327781.84, 246168.69]
    rvir = 2177.5625

    #GU = GadgetPhysical(read_header(snap_base))

    read_particles_in_volume(snap_base, blocks, gpos, 2rvir, 
                            use_keys=use_keys)#, GU
end
# Load the data
sim_path = "/gpfs/work/pn68va/di67meg/LocalUniverse/"
snap = 36

data_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/radial_B/box/"

function run()

    folders = ["zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_inj"]
    data_paths = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/radial_B/" .* folders .* "/"

    sim_paths = "/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/" .* [
        "mhd_cr8p24eDpp_1e-17/", "mhd_cr8p24eDpp_5e-17/", "mhd_cr8p24e/"]

    snaps = [74, 12, 12]

    for j = 1:length(data_paths)

        data = read_data(sim_paths[j], snaps[j], false)

        Bfield_functions = [Bfield_sim, Bfield_Beta, Bfield_vturb, Bfield_FF, Bfield_dyn_l, Bfield_dyn_h]
        Bfield_names = ["B_sim", "B_beta", "B_vturb", "B_FF", "B_dyn_l", "B_dyn_h"]

        for i = 1:length(Bfield_functions)
            write_binning(data_paths[j] * Bfield_names[i] * ".dat", data, Bfield_functions[i])
        end

    end

end

run()

folders = [#"zoom_dpp_1e-17", "zoom_dpp_5e-17", 
            "zoom_inj"]
data_paths = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/radial_B/" .* folders .* "/"

sim_paths = "/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/" .* [
    #"mhd_cr8p24eDpp_1e-17/", "mhd_cr8p24eDpp_5e-17/", 
    "mhd_cr8p24e/"]

snaps = [#74, 12, 
        12]

const global GU = GadgetPhysical(read_header("/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/mhd_cr8p24eDpp_1e-17/snapdir_074/snap_074.0"))

for j = 1:length(data_paths)

    data = read_data(sim_paths[j], snaps[j], false)

    Bfield_functions = [Bfield_sim, Bfield_Beta, Bfield_vturb, Bfield_FF, Bfield_dyn_l, Bfield_dyn_h]
    Bfield_names = ["B_sim", "B_beta", "B_vturb", "B_FF", "B_dyn_l", "B_dyn_h"]

    for i = 1:length(Bfield_functions)
        write_binning(data_paths[j] * Bfield_names[i] * ".dat", data, Bfield_functions[i])
    end

end