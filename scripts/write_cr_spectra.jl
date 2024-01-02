using GadgetIO, GadgetUnits
using Printf
using SpectralCRsUtility
using PyPlot, PyPlotUtility
using ProgressMeter
using Base.Threads
using DelimitedFiles


struct CR
    spec::CRMomentumDistribution
    j_nu::Vector{Float64}
end


function get_synch_power(data, h, i;
    pmin=1.e-1, pmax=1.e5)

    GU = GadgetPhysical(h)

    # define observation frequency in Hz and shift to redshift
    Nnu = 50

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)

    ν_arr = 10.0 .^ (LinRange(7, 10, Nnu))
    j_ν = Vector{Float64}(undef, Nnu)

    for iNu = 1:Nnu

        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])

        B = √(data["BFLD"][1, i]^2 +
              data["BFLD"][2, i]^2 +
              data["BFLD"][3, i]^2)

        j_ν[iNu] = synchrotron_emission(norm, slope, cut, B, par, ν0=ν_arr[iNu],
            reduce_spectrum=true,
            integrate_pitch_angle=true)

    end # Nu

    return j_ν
end


function read_data(snap)

    sim_path = "/mnt/home/lboess/ceph/LocalUniverse/Coma/L5/mhd_cr8p24eDpp/"
    snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
    save_file = sim_path * "read_positions/max_CReP_ids_$(@sprintf("%03i", snap)).dat"

    read_positions = load_read_positions(save_file)
    blocks = ["BFLD", "CReN", "CReS", "CReC"]

    return read_blocks_filtered(snap_base, blocks; read_positions), read_header(snap_base)
end


function construct_spectrum(data, i;
    pmin=1.e-1, pmax=1.e5)

    return CRMomentumDistribution(data["CReN"][:, i], data["CReS"][:, i], data["CReC"][i],
        pmin, pmax, 1.0, 3)
end


function convert_data(data, h, i)
    return CR(construct_spectrum(data, i), get_synch_power(data, h, i))
end




# function write_spectra(spectra, id)
#     Nfiles = size(spectra, 2)
#     filename = "/mnt/home/lboess/ceph/LocalUniversePapers/zooms/coma/data/spectra_$id.dat"
#     f = open(filename)
#     write(f, Nfiles)
#     for Nfile ∈ 1:Nfiles
#         write(f, spectra[Nfile].spec.bound)
#         write(f, spectra[Nfile].spec.norm)
#         write(f, spectra[Nfile].j_nu)
#     end
#     close(f)
# end

# function run(snap_range)

#     ids = [ 549847351403, 549834206059, 549816871198, 549822811914, 549834873392, 549850851667,
#         549837820882, 549784127556, 549848445311, 549782540344, 549813136996, 549890489504,
#         549836262508, 549849186219, 549888791935, 549815238762, 549849858007, 549771151247,
#         549778132736, 549839258204, 549835496070, 549856547328, 549839327739 ]

#     for (Nid, id) ∈ enumerate(ids)

#         spectra = Vector{CR}(undef, length(snap_range))

#         for (i, snap) ∈ enumerate(snap_range)

#             @info "snap $snap"
#             data, header = read_data(snap)
#             spectra[i] = convert_data(data, header, Nid)
#             data = nothing
#             header = nothing
#         end

#         write_spectra(spectra, id)
#     end
# end

function write_spectra(spectra, ids)

    Nids = length(ids)
    Nfiles = size(spectra, 2)
    filename = "/mnt/home/lboess/ceph/LocalUniversePapers/zooms/coma/data/spectra.dat"
    f = open(filename, "w")
    write(f, Nids)
    write(f, ids)
    write(f, Nfiles)
    for Nfile ∈ 1:Nfiles, Nid ∈ 1:Nids
        write(f, spectra[Nid, Nfile].spec.bound)
        write(f, spectra[Nid, Nfile].spec.norm)
        write(f, spectra[Nid, Nfile].j_nu)
    end
    close(f)

end

function run(snap_range)

    @info "running on $(nthreads()) threads"

    ids = [549847351403, 549834206059, 549816871198, 549822811914, 549834873392, 549850851667,
        549837820882, 549784127556, 549848445311, 549782540344, 549813136996, 549890489504,
        549836262508, 549849186219, 549888791935, 549815238762, 549849858007, 549771151247,
        549778132736, 549839258204, 549835496070, 549856547328, 549839327739]


    spectra = Matrix{CR}(undef, 23, length(snap_range))

    for (i, snap) ∈ enumerate(snap_range)
        @info "snap $snap"
        data, header = read_data(snap)

        @threads for Nid ∈ 1:length(ids)
            spectra[Nid, i] = convert_data(data, header, Nid)
        end
    end

    write_spectra(spectra, ids)
end

#run(0:73)

function write_times(snap_range)

    t = Vector{Float64}(undef, length(snap_range))

    for (i, snap) ∈ enumerate(snap_range)
        @info "snap $snap"
        sim_path = "/mnt/home/lboess/ceph/LocalUniverse/Coma/L5/mhd_cr8p24eDpp/"
        snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"

        h = read_header(snap_base)
        t[i] = age(h)
    end

    fo = "/mnt/home/lboess/ceph/LocalUniversePapers/zooms/coma/data/times.txt"
    writedlm(fo, t)
end

write_times(0:73)