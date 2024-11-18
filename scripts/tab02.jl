include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "shared.jl"))

using Printf

phase_map_path = data_path * "phase_maps/"


for folder in ["box", "zoom_inj", "zoom_dpp"]
    for B in ["B_sim", "B_beta50", "B_Pturb", "B_FF", "B_dyn_l", "B_dyn_h"]
        filename = phase_map_path * "$folder/bin_1D_synch_emissivity_144MHz_$B.dat"
        ne, syn = read_1D_data(filename)
        println(folder, " ", B, ": ", syn[46])
    end
end
