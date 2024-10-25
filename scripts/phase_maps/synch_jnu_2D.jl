println("allocating cores")
using Distributed, ClusterManagers
# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        N_tasts_ = 4
        println("allocating $N_tasts_ normal tasks")
        addprocs(N_tasts_)
    end
end

println("done")
flush(stdout);
flush(stderr);

println("loading packages")
@everywhere using GadgetIO, GadgetUnits
@everywhere using Printf
@everywhere using SpectralCRsUtility
@everywhere using Base.Threads
@everywhere using Statistics
println("done")
flush(stdout);
flush(stderr);

# if length(ARGS) == 1
#     @everywhere global const Bfield_flag = 1
#     @everywhere global const flux = false
# else
if ARGS[1] == "sim"
    @everywhere const global Bfield_flag = 1
elseif ARGS[1] == "beta"
    @everywhere const global Bfield_flag = 2
elseif ARGS[1] == "vturb"
    @everywhere const global Bfield_flag = 3
elseif ARGS[1] == "ff"
    @everywhere const global Bfield_flag = 4
elseif ARGS[1] == "dyn_l"
    @everywhere const global Bfield_flag = 5
elseif ARGS[1] == "dyn_h"
    @everywhere const global Bfield_flag = 6
end

# if ARGS[2] == "flux"
#     @everywhere global const flux = true 
# else
@everywhere const global flux = false
# end
#end

# mapping settings
# @everywhere const snap = 36
# @everywhere const global sim_path = "path/to/simulation/"
# @everywhere const data_path = "path/to/data"
# @everywhere const snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
@everywhere const global snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"

@everywhere const map_path = "/e/ocean2/users/lboess/PaperRepos/2310.13734/data/phase_maps/box/"

@everywhere const global GU = GadgetPhysical(GadgetIO.read_header(snap_base))

@everywhere const global center_comov = [247.980, 245.480, 255.290] .* 1.e3
@everywhere const global center = center_comov .* GU.x_physical
@everywhere const global radius_limits = [10_000.0, 240_000.0 * GU.x_physical] # for flux

@everywhere include(joinpath(@__DIR__, "bin_2D.jl"))
@everywhere include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))
@everywhere include(joinpath(@__DIR__, "..", "allsky", "synchrotron.jl"))

const jobs = RemoteChannel(() -> Channel{Int}(32))
const results = RemoteChannel(() -> Channel{Matrix}(80000))

# phase map settings
@everywhere const x_lim = [1.e-9, 1.0]
@everywhere const y_lim = [1.e-52, 1.e-37]

@everywhere const Nbins = 100

@everywhere function get_phase_map(rho, j_nu)
    bin_2D_log(rho, j_nu, x_lim, y_lim,
        Nbins=Nbins, show_progress=false, calc_mean=false)
end

@everywhere function synch_phase_map_of_subfile(subfile, Bfield_function, Btype)

    println("syn (B$Btype): subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["POS", "MASS", "RHO", "U", "BFLD", "VRMS",
        "CReN", "CReS", "CReC"])

    #println("\ttemperature and density")
    flush(stdout)
    flush(stderr)
    ne = data["RHO"] .* GU.rho_ncm3

    # calculate synchrotron emissivity
    j_nu = get_synchrotron(data, 144.e6, Bfield_function, flux)

    j_nu[iszero.(j_nu)] .= 1.e-100
    j_nu[j_nu.<1.e-100] .= 1.e-100
    j_nu[isnan.(j_nu)] .= 1.e-100
    j_nu[isinf.(j_nu)] .= 1.e-100

    # calculate synchrotron power
    unit = "erg/cm^3/s/Hz"

    println("\tsynch done!\n\tminimum = $(minimum(j_nu[j_nu .> 0.0])) $unit\n\tmaximum = $(maximum(j_nu)) $unit\n\tsum = $(sum(j_nu)) $unit")
    flush(stdout)
    flush(stderr)

    phase_P_nu = get_phase_map(ne, j_nu)
    j_nu = nothing
    data = nothing
    ne = nothing
    GC.gc()

    return phase_P_nu
end

@everywhere function synch_phase_map_of_subfile(subfile)
    if Bfield_flag == 1
        return synch_phase_map_of_subfile(subfile, Bfield_sim, "sim")
    elseif Bfield_flag == 2
        return synch_phase_map_of_subfile(subfile, Bfield_Beta, "beta")
    elseif Bfield_flag == 3
        return synch_phase_map_of_subfile(subfile, Bfield_vturb, "vturb")
    elseif Bfield_flag == 4
        return synch_phase_map_of_subfile(subfile, Bfield_FF, "ff")
    elseif Bfield_flag == 5
        return synch_phase_map_of_subfile(subfile, Bfield_dyn_l, "dyn_l")
    elseif Bfield_flag == 6
        return synch_phase_map_of_subfile(subfile, Bfield_dyn_h, "dyn_h")
    end
end

@everywhere function get_synch_phase_map_filename()

    if Bfield_flag == 1
        return map_path * "phase_map_synch_emissivity_144MHz_B_sim.dat"
    elseif Bfield_flag == 2
        return map_path * "phase_map_synch_emissivity_144MHz_B_beta50.dat"
    elseif Bfield_flag == 3
        return map_path * "phase_map_synch_emissivity_144MHz_B_Pturb.dat"
    elseif Bfield_flag == 4
        return map_path * "phase_map_synch_emissivity_144MHz_B_FF.dat"
    elseif Bfield_flag == 5
        return map_path * "phase_map_synch_emissivity_144MHz_B_dyn_l.dat"
    elseif Bfield_flag == 6
        return map_path * "phase_map_synch_emissivity_144MHz_B_dyn_h.dat"
    end

end


@everywhere function write_phase_map(filename, phase_map)
    Nbins = size(phase_map, 1)
    f = open(filename, "w")
    write(f, Nbins)
    write(f, x_lim)
    write(f, y_lim)
    write(f, phase_map)
    close(f)
end


@everywhere function do_work(jobs, results)
    while true
        job_id = take!(jobs)
        put!(results, synch_phase_map_of_subfile(job_id - 1))
    end
end


function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end


function run_phase_maps(filename)

    println("starting workers")

    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results)
    end

    n = 2048
    #n = 4

    sum_phase_P_nu = zeros(Float64, Nbins, Nbins)

    errormonitor(@async make_jobs(n)) # feed the jobs channel with "n" jobs

    println("running")
    flush(stdout)
    flush(stderr)

    @time while n > 0 # print out results

        phase_P_nu = take!(results)

        # sum up contribution
        sum_phase_P_nu += phase_P_nu
        n -= 1

        phase_P_nu = nothing
    end
    flush(stdout)
    flush(stderr)

    println("Pnu phase map")
    println("    Min:    $(minimum(sum_phase_P_nu))")
    println("    Max:    $(maximum(sum_phase_P_nu))")
    println("    Mean:   $(mean(sum_phase_P_nu))")
    println("    Median: $(median(sum_phase_P_nu))")

    flush(stdout)
    flush(stderr)
    write_phase_map(filename, sum_phase_P_nu)
end

filename = get_synch_phase_map_filename()
run_phase_maps(filename)
