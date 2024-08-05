using GadgetIO, GadgetUnits 
using SPHKernels, SPHtoGrid
using Printf
using Base.Threads
using LinearAlgebra

const global GU = GadgetPhysical(read_header("/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/mhd_cr8p24eDpp_1e-17/snapdir_074/snap_074.0"))

include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))

const global gpos = [244977.58, 327853.62, 245989.75] .* 1.e-3
const global side_length = 5_000.0 

function read_data(snap_base, use_keys)

    blocks = ["POS", "VEL", "HSML", "RHO", "U", "MASS", "VRMS", "BFLD"]

    # coma
    cube = GadgetCube(gpos .* 1.e3, side_length)

    return read_particles_in_geometry(snap_base, blocks, cube, 
                    use_keys=use_keys)
end

function rescale_bfield(data, Bfield_funtion)
    
    B = Matrix{Float32}(undef, 3, size(data["BFLD"], 2))
    
    @threads for i ∈ 1:size(data["BFLD"], 2)
        Bnorm = norm(data["BFLD"][:, i])
        B_ref = Bfield_funtion(data, i)
        
        if iszero(Bnorm)
            B_rescale = 0.0
        else
            B_rescale = B_ref / Bnorm
        end

        for j ∈ 1:3
            B[j, i] = data["BFLD"][j, i] * Float32(B_rescale)
        end
    end

    return B
end

function write_snapshot(data, B)

    Npart = size(data["POS"], 2)

    blocks = ["POS", "VEL", "HSML", "RHO", "U", "MASS"]
    h = read_header("/gpfs/work/pn68va/di67meg/LocalUniverse/snapdir_036/snap_036.0")
    h.npart[1] = Npart 
    h.npart[2:end] .= 0
    h.nall = UInt32.(h.npart)
    h.npartTotalHighWord .= 0
    h.num_files = 1
    h.massarr .= 0.0

    f = open("snap_000", "w")
    write_header(f, h)
    for block ∈ blocks 
        write_block(f, data[block], block)
    end
    write_block(f, B, "BFLD")
    close(f)

end


function run_it(snap_base, Bfield_funtion, map_filename, use_keys)

    println("reading data")
    data = read_data(snap_base, use_keys)
    println("rescaling B")
    B = rescale_bfield(data, Bfield_funtion)
    println("writing snapshot")
    write_snapshot(data, B)
    println("running smac")

    h = read_header("snap_000")
    write_smac1_par(SNAP_FILE="snap_000", KERNEL_TYPE=2,
        IMG_XY_SIZE = side_length, IMG_Z_SIZE=side_length, IMG_SIZE=1024,
        OUTPUT_MAP=12, OUTPUT_SUB=0,
        PROJECT=2,
        CENTER_X=gpos[1], CENTER_Y=gpos[2], CENTER_Z=gpos[3],
        PREFIX_OUT=map_filename, h=h)

        # execute smac
    run(`./Smac_6.1 smac1.par`)
end


Bfield_functions = [Bfield_sim, Bfield_Beta, 
    Bfield_vturb,
    Bfield_FF, Bfield_dyn_l, Bfield_dyn_h]
Bfield_names = ["B_sim", "B_beta", 
    "B_vturb_rescaled",
    "B_FF_rescaled", "B_dyn_l_rescaled", "B_dyn_h_rescaled"]

# box 


folders = ["zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_inj"]
data_paths = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/coma/" .* folders .* "/"

sim_paths = "/gpfs/work/pn68va/di67meg/LocalUniverseZoom/Coma/L5/" .* [
    "mhd_cr8p24eDpp_1e-17/", "mhd_cr8p24eDpp_5e-17/", "mhd_cr8p24e/"]

snaps = [74, 12, 12]

for j = 1:1#length(data_paths)

    for i = 5:5#4:length(Bfield_functions)
        map_filename = "coma_5Mpc_$(@sprintf("%03i", snaps[j])).RM_" * Bfield_names[i]
        println(map_filename)
        run_it(sim_paths[j] * "snapdir_$(@sprintf("%03i", snaps[j]))/snap_$(@sprintf("%03i", snaps[j]))",
            Bfield_functions[i], map_filename, false)
    end

end
