println("allocating cores")
using Distributed, ClusterManagers
# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        N_tasts_ = 2
        println("allocating $N_tasts_ normal tasks")
        addprocs(N_tasts_)
    end
end

println("done")
flush(stdout); flush(stderr)

println("loading packages")
@everywhere using GadgetIO, GadgetUnits
@everywhere using Printf
@everywhere using SpectralCRsUtility
@everywhere using Base.Threads
@everywhere using Statistics
println("done")
flush(stdout); flush(stderr)


if ARGS[1] == "sim"
    @everywhere global const Bfield_flag = 1
elseif ARGS[1] == "beta"
    @everywhere global const Bfield_flag = 2
elseif ARGS[1] == "vturb"
    @everywhere global const Bfield_flag = 3
elseif ARGS[1] == "ff"
    @everywhere global const Bfield_flag = 4
elseif ARGS[1] == "dyn_l"
    @everywhere global const Bfield_flag = 5
elseif ARGS[1] == "dyn_h"
    @everywhere global const Bfield_flag = 6
end


# mapping settings
@everywhere const snap = 12
# @everywhere const global sim_path = "/e/ocean3/Local/3072/nonrad_mhd_crs/"
@everywhere const global sim_path = "/e/ocean2/users/lboess/LocalUniverseZooms/L5/mhd_cr8p24eDpp_5e-17/"
@everywhere const snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
@everywhere const map_path = joinpath(@__DIR__, "..", "..", "data", "phase_maps", "zoom_dpp_high")

@everywhere global const GU = GadgetPhysical(GadgetIO.read_header(snap_base))

@everywhere include(joinpath(@__DIR__, "bin_2D.jl"))
@everywhere include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))
@everywhere include(joinpath(@__DIR__, "..", "allsky", "synchrotron.jl"))

const jobs = RemoteChannel(()->Channel{Int}(32))
const results = RemoteChannel(()->Channel{Matrix}(80000))

# phase map settings
@everywhere const x_lim = [1.e-9, 1.0]
@everywhere const y_lim = [10.0,  1.e10]

@everywhere const Nbins = 100

@everywhere function get_phase_map(rho, T, quantity)
    bin_2D_log( rho, T, x_lim, y_lim, quantity, 
                Nbins = Nbins, show_progress = false, calc_mean=false)
end

@everywhere function synch_phase_map_of_subfile(subfile, Bfield_function, Btype)

    println("syn (B$Btype): subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["POS", "MASS", "RHO", "U", "BFLD", "VRMS", 
                             "CReN", "CReS", "CReC"])

    println("\ttemperature and density")
    flush(stdout); flush(stderr)
    T  = data["U"]   .* GU.T_K 
    ne = data["RHO"] .* GU.rho_ncm3

    # calculate synchrotron emissivity
    j_nu = get_synchrotron(data, 144.e6, Bfield_function)

    # calculate synchrotron power
    P_nu = @. j_nu * (data["MASS"] .* GU.m_cgs) / (data["RHO"] .* GU.rho_cgs) * 1.e-7

    println("\tsynch done!\n\tmaximum = $(maximum(P_nu)) W/Hz\n\tsum = $(sum(P_nu)) W/Hz")
    flush(stdout)
    flush(stderr)

    phase_P_nu = get_phase_map(ne, T, P_nu)
    j_nu = nothing
    P_nu = nothing
    data = nothing
    T = ne = nothing
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
        return map_path * "/bin_2D_synch_power_144MHz_B_sim.dat"
    elseif Bfield_flag == 2
        return map_path * "/bin_2D_synch_power_144MHz_B_beta50.dat"
    elseif Bfield_flag == 3
        return map_path * "/bin_2D_synch_power_144MHz_B_01Pturb.dat"
    elseif Bfield_flag == 4
        return map_path * "/bin_2D_synch_power_144MHz_B_FF.dat"
    elseif Bfield_flag == 5
        return map_path * "/bin_2D_synch_power_144MHz_B_dyn_l.dat"
    elseif Bfield_flag == 6
        return map_path * "/bin_2D_synch_power_144MHz_B_dyn_h.dat"
    end

end

@everywhere function write_phase_map(filename, phase_map)
    Nbins = size(phase_map,1)
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
        put!(results, synch_phase_map_of_subfile(job_id-1))
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
    
    n = 16
    #n = 4

    sum_phase_P_nu = zeros(Float64, Nbins, Nbins) 

    errormonitor(@async make_jobs(n)); # feed the jobs channel with "n" jobs
    
    println("running")
    flush(stdout); flush(stderr)

    @time while n > 0 # print out results

        phase_P_nu = take!(results)

        # sum up contribution
        sum_phase_P_nu += phase_P_nu
        n -= 1

        phase_P_nu = nothing
    end
    flush(stdout); flush(stderr)

    println("Pnu phase map")
    println("    Min:    $(minimum(sum_phase_P_nu))")
    println("    Max:    $(maximum(sum_phase_P_nu))")
    println("    Mean:   $(mean(sum_phase_P_nu))")
    println("    Median: $(median(sum_phase_P_nu))")

    flush(stdout); flush(stderr)
    write_phase_map(filename, sum_phase_P_nu)
end

filename = get_synch_phase_map_filename()
run_phase_maps(filename)
