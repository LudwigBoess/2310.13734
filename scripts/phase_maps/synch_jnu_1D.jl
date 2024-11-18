println("allocating cores")
using Distributed, ClusterManagers

# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks, using $(2 * parse(Int64, ENV["SLURM_CPUS_PER_TASK"])) threads each")
    withenv("JULIA_NUM_THREADS" => 2 * parse(Int64, ENV["SLURM_CPUS_PER_TASK"])) do
        addprocs(addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))) # spawn 4 workers with 2 threads each
    end
catch err
    if isa(err, KeyError)
        println("allocating 4 normal tasks")
        addprocs(4)
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
@everywhere using ProgressMeter
# @everywhere using ThreadPinning
# pinthreads(:cores)

println("done")
flush(stdout);
flush(stderr);

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


#@everywhere const snap = 74
#@everywhere const global sim_path = "/gpfs/work/pn36ze/di93son/LocalUniverse/Coma/L5/cr6p20eDpp/"
#@everywhere const snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
@everywhere const global snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"
#@everywhere const global snap_base = "/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20e/snapdir_074/snap_074"
@everywhere const global GU = GadgetPhysical(GadgetIO.read_header(snap_base))

@everywhere const map_path = joinpath(@__DIR__, "..", "..", "data", "phase_maps", "box")

@everywhere global const center_comov = zeros(3)
@everywhere const global radius_limits = [0.0, Inf]

@everywhere include(joinpath(@__DIR__, "bin_1D.jl"))
@everywhere include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))
@everywhere include(joinpath(@__DIR__, "..", "allsky", "synchrotron.jl"))

const jobs = RemoteChannel(() -> Channel{Int}(32))
const results = RemoteChannel(() -> Channel{Tuple}(1600))

# phase map settings
@everywhere const x_lim = [1.e-9, 1.0]
@everywhere const Nbins = 100


@everywhere function synch_1D_bin_of_subfile(subfile, Bfield_function, Btype)

    println("syn (B$Btype): subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["POS", "MASS", "RHO", "U", "BFLD", "VRMS",
                             "CReN", "CReS", "CReC"])

    println("\ttemperature and density")
    flush(stdout)
    flush(stderr)
    ne = data["RHO"] .* GU.rho_ncm3

    # calculate synchrotron emissivity
    j_nu = get_synchrotron(data, 144.e6, Bfield_function, false)

    j_nu[iszero.(j_nu)] .= 1.e-100
    j_nu[j_nu.<1.e-100] .= 1.e-100
    j_nu[isnan.(j_nu)] .= 1.e-100
    j_nu[isinf.(j_nu)] .= 1.e-100

    # only bin the values within the phase range
    sel = findall(j_nu .> 1.e-52)
    if length(sel) > 0
        j_nu = j_nu[sel]
        ne = ne[sel]

        unit = "erg/cm^3/s/Hz"

        println("\tsynch done!\n\tminimum = $(minimum(j_nu[j_nu .> 0.0])) $unit\n\tmaximum = $(maximum(j_nu)) $unit\n\tsum = $(sum(j_nu)) $unit")
        flush(stdout)
        flush(stderr)

        ne_count, jnu_sum = bin_1D_log(ne, x_lim, j_nu, show_progress=false, calc_sigma=false)
    else
        ne_count = zeros(Int64, Nbins)
        jnu_sum = zeros(Float64, Nbins)
    end
    j_nu = nothing
    data = nothing
    ne = nothing
    sel = nothing
    GC.gc()

    return ne_count, jnu_sum
end



@everywhere function write_synch_binning(filename, jnu, jnu_sigma)
    Nbins = length(jnu)
    f = open(filename, "w")
    write(f, Nbins)
    write(f, x_lim)
    write(f, jnu)
    write(f, jnu_sigma)
    close(f)
end


@everywhere function synch_1D_bin_of_subfile(subfile)
    if Bfield_flag == 1
        return synch_1D_bin_of_subfile(subfile, Bfield_sim, "sim")
    elseif Bfield_flag == 2
        return synch_1D_bin_of_subfile(subfile, Bfield_Beta, "beta")
    elseif Bfield_flag == 3
        return synch_1D_bin_of_subfile(subfile, Bfield_vturb, "vturb")
    elseif Bfield_flag == 4
        return synch_1D_bin_of_subfile(subfile, Bfield_FF, "ff")
    elseif Bfield_flag == 5
        return synch_1D_bin_of_subfile(subfile, Bfield_dyn_l, "dyn_l")
    elseif Bfield_flag == 6
        return synch_1D_bin_of_subfile(subfile, Bfield_dyn_h, "dyn_h")
    end
end


@everywhere function get_synch_phase_map_filename()

    if Bfield_flag == 1
        return map_path * "/bin_1D_synch_emissivity_144MHz_B_sim.dat"
    elseif Bfield_flag == 2
        return map_path * "/bin_1D_synch_emissivity_144MHz_B_beta50.dat"
    elseif Bfield_flag == 3
        return map_path * "/bin_1D_synch_emissivity_144MHz_B_Pturb.dat"
    elseif Bfield_flag == 4
        return map_path * "/bin_1D_synch_emissivity_144MHz_B_FF.dat"
    elseif Bfield_flag == 5
        return map_path * "/bin_1D_synch_emissivity_144MHz_B_dyn_l.dat"
    elseif Bfield_flag == 6
        return map_path * "/bin_1D_synch_emissivity_144MHz_B_dyn_h.dat"
    end

end


@everywhere function do_work(jobs, results) # define work function everywhere
    while true
        job_id = take!(jobs)
        put!(results, synch_1D_bin_of_subfile(job_id - 1))
    end
end

function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end


function run_jnu_1D()

    println("starting workers")

    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results)
    end

    n = 2048
    #n = 4

    sum_ne_count = zeros(Int64, Nbins)
    sum_jnu = zeros(Float64, Nbins)

    errormonitor(@async make_jobs(n)) # feed the jobs channel with "n" jobs

    println("running")
    flush(stdout)
    flush(stderr)

    @time while n > 0 # print out results

        ne_count, jnu = take!(results)
        sum_ne_count .+= ne_count
        sum_jnu .+= jnu
        n -= 1

        ne_count = jnu = nothing
        GC.gc()
    end
    flush(stdout)
    flush(stderr)

    jnu_sigma = σ_1D_quantity(sum_jnu, sum_ne_count)

    for i = 1:length(sum_ne_count)
        if !iszero(sum_ne_count[i])
            sum_jnu[i] /= sum_ne_count[i]
        end
    end

    println("jnu 1D")
    println("    Min:    $(minimum(sum_jnu))")
    println("    Max:    $(maximum(sum_jnu))")
    println("    Mean:   $(mean(sum_jnu))")
    println("    Median: $(median(sum_jnu))")

    flush(stdout)
    flush(stderr)
    filename = get_synch_phase_map_filename()
    write_synch_binning(filename, sum_jnu, jnu_sigma)
end

run_jnu_1D()