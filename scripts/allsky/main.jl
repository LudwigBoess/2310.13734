println("allocating cores")
using Distributed, ClusterManagers

# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks, using $(2 * parse(Int64, ENV["SLURM_CPUS_PER_TASK"])) threads each")
    withenv("JULIA_NUM_THREADS" => 2 * parse(Int64, ENV["SLURM_CPUS_PER_TASK"])) do
        addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"])) # spawn 4 workers with 2 threads each
    end
catch err
    if isa(err, KeyError)
        println("allocating 20 normal tasks")
        addprocs(20)
    end
end

println("pinning threads")
@everywhere using ThreadPinning
distributed_pinthreads(:cores)
println(distributed_getcpuids())
println(distributed_getispinned())
println("done")
flush(stdout);
flush(stderr);

println("loading packages")
@everywhere using GadgetIO, GadgetUnits
@everywhere using Printf
@everywhere using SPHKernels, SPHtoGrid
@everywhere using SpectralCRsUtility
@everywhere using Base.Threads
@everywhere using Statistics
@everywhere using Healpix

println("done")
flush(stdout);
flush(stderr);

if length(ARGS) == 1
    @everywhere const global Bfield_flag = 1
else

    if ARGS[2] == "sim"
        @everywhere const global Bfield_flag = 1
    elseif ARGS[2] == "beta"
        @everywhere const global Bfield_flag = 2
    elseif ARGS[2] == "vturb"
        @everywhere const global Bfield_flag = 3
    elseif ARGS[2] == "ff"
        @everywhere const global Bfield_flag = 4
    elseif ARGS[2] == "dyn_l"
        @everywhere const global Bfield_flag = 5
    elseif ARGS[2] == "dyn_h"
        @everywhere const global Bfield_flag = 6
    end

end

# include map functions
@everywhere include("config.jl")
@everywhere include("CReE.jl")
@everywhere include("Bfld.jl")
@everywhere include("synchrotron.jl")

println("Bfield_flag = ", Bfield_flag)

"""
    run_it()

Main function
"""
function run_it()


    if ARGS[1] == "CReE_seed"

        filename = map_path * "allsky_CReP_seed_$viewpoint.fits"
        # `distributed_allsky_map` is imported from SPHtoGrid.jl
        distributed_allsky_map(filename, Nside, Nfiles, CReE_seed_maps_of_subfile)

    elseif ARGS[1] == "CReE"

        filename = map_path * "allsky_CReE_5Mpc_$viewpoint.fits"
        distributed_allsky_map(filename, Nside, Nfiles, CReE_maps_of_subfile, reduce_image=false)

    elseif ARGS[1] == "B"

        filename = get_B_filename()
        distributed_allsky_map(filename, Nside, Nfiles, Bfld_map_of_subfile)

    elseif ARGS[1] == "syn"

        filename = get_synch_filename()
        distributed_allsky_map(filename, Nside, Nfiles, synch_maps_of_subfile, reduce_image=false)

    end

end

run_it()