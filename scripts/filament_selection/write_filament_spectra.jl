using GadgetIO, GadgetUnits
using SpectralCRsUtility
using Printf
using ProgressMeter
using DelimitedFiles

include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))

const global h = read_header("/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20e/snapdir_074/snap_074.0")
const global GU = GadgetPhysical(h)

function read_data(snap_base, use_keys)

    blocks = ["POS", "MASS", "HSML", "RHO", "U", 
              "BFLD", "VRMS", "CReN", "CReS", "CReC"]

    # found by selecting particles by hand
    cylinder = GadgetCylinder([247.75, 332.0, 248.5] .* 1.e3, 
                              [248.5, 332.75, 249.25] .* 1.e3, 
                              2.5e3)

    return read_particles_in_geometry(snap_base, blocks, cylinder; use_keys)
end

"""
    get_synchrotron(data, nu, Bfield_function, show_progress=false)

Calculate synchrotron emissivity for a given data set at observational frequency `nu` for a magnetic field defined by `Bfield_function`.
"""
function get_synchrotron(data, nu, Bfield_function)

    Npart = length(data["CReC"])
    j_ν = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(1.0, 1.e5, Nbins)

    @threads for i ∈ eachindex(j_ν)

        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]

        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])

        B = Bfield_function(data, i)

        j_ν[i] = synchrotron_emission(norm, slope, cut, B, par, ν0=nu,
            reduce_spectrum=true,
            integrate_pitch_angle=true,
            convert_to_mJy=false)

        m = data["MASS"][i] * GU.m_cgs
        rho = data["RHO"][i] * GU.rho_cgs

        # convert to [W/Hz]
        j_ν[i] *= m / rho * 1.e-7
    end

    sum(j_ν)
end


function get_synch_spectrum(data, Bfield_function)

    @info "setting up spectra"
    # define observation frequency in Hz
    Nnu = 50
    ν_arr = 10.0 .^ (LinRange(7, 10, Nnu))
    P = Vector{Float64}(undef, Nnu)

    # cr setup 
    @showprogress for iNu = 1:Nnu
        # get total synch power per frequency
        P[iNu] = get_synchrotron(data, ν_arr[iNu], Bfield_function)
    end # Nu

    @info "done"
    flush(stdout)
    flush(stderr)

    return ν_arr, P
end

function get_energy_spectrum(data)

    Nbins = size(data["CReN"], 1)
    Npart = length(data["CReC"])
    N = Vector{Float64}(undef, Nbins)

    par = CRMomentumDistributionConfig(1.0, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    @showprogress for i ∈ 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])
        rho = data["RHO"][i] * GU.rho_cgs

        @threads for j = 1:Nbins
            if bounds[j] > cut
                continue
            end
            hbound_proper = bounds[j+1] > cut ? cut : bounds[j+1]
            N[j] += SpectralCRsUtility.density_integral(bounds[j], hbound_proper, norm[j], slope[j], rho)
        end

    end

    centers = momentum_bin_centers(par)

    return centers, N
end

function get_e_spectrum(data)

    Nbins = size(data["CReN"], 1)
    Npart = length(data["CReC"])
    E = Vector{Float64}(undef, Nbins)

    par = CRMomentumDistributionConfig(1.0, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    @showprogress for i ∈ 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])
        rho = data["RHO"][i] * GU.rho_cgs

        @threads for j = 1:Nbins
            if bounds[j] > cut
                continue
            end
            hbound_proper = bounds[j+1] > cut ? cut : bounds[j+1]
            E[j] += SpectralCRsUtility.energy_integral(bounds[j], hbound_proper, norm[j], slope[j], rho)
        end

    end

    centers = momentum_bin_centers(par)

    return centers, E
end

function get_f_spectrum(data)

    Nbins = size(data["CReN"], 1)
    Npart = length(data["CReC"])
    f = Vector{Float64}(undef, Nbins)

    par = CRMomentumDistributionConfig(1.0, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    @showprogress for i ∈ 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])
        rho = data["RHO"][i] * GU.rho_cgs

        @threads for j = 1:Nbins
            if bounds[j] > cut
                continue
            end
            #hbound_proper = bounds[j+1] > cut ? cut : bounds[j+1]
            f[j] += norm[j] #SpectralCRsUtility.energy_integral(bounds[j], hbound_proper, norm[j], slope[j], rho)
        end

    end

    centers = momentum_bin_centers(par)

    return centers, f
end

function run(snap_base, use_keys, out_dir)

    @info "reading data"
    data = read_data(snap_base, use_keys)

    # @info "calculating energy spectrum"
    # centers, E = get_energy_spectrum(data)
    # fo = joinpath(out_dir, "energy_spectrum.dat")
    # writedlm(fo, [centers E])

    @info "calculating energy spectrum"
    centers, E = get_f_spectrum(data)
    fo = joinpath(out_dir, "f_spectrum.dat")
    writedlm(fo, [centers E])

    @info "calculating synchrotron spectrum"

    # @info "Bfield: sim"
    # ν_arr, P = get_synch_spectrum(data, Bfield_sim)
    # fo = joinpath(out_dir, "synch_spectrum_sim.dat")
    # writedlm(fo, [ν_arr P])

    # @info "Bfield: ff"
    # ν_arr, P = get_synch_spectrum(data, Bfield_FF)
    # fo = joinpath(out_dir, "synch_spectrum_ff.dat")
    # writedlm(fo, [ν_arr P])

    # @info "Bfield: beta"
    # ν_arr, P = get_synch_spectrum(data, Bfield_Beta)
    # fo = joinpath(out_dir, "synch_spectrum_beta.dat")
    # writedlm(fo, [ν_arr P])

    @info "Bfield: Bfield_vturb"
    ν_arr, P = get_synch_spectrum(data, Bfield_vturb)
    fo = joinpath(out_dir, "synch_spectrum_vturb1.dat")
    writedlm(fo, [ν_arr P])

    # @info "Bfield: Bfield_dyn_l"
    # ν_arr, P = get_synch_spectrum(data, Bfield_dyn_l)
    # fo = joinpath(out_dir, "synch_spectrum_dyn_l.dat")
    # writedlm(fo, [ν_arr P])

    # @info "Bfield: Bfield_dyn_h"
    # ν_arr, P = get_synch_spectrum(data, Bfield_dyn_h)
    # fo = joinpath(out_dir, "synch_spectrum_dyn_h.dat")
    # writedlm(fo, [ν_arr P])
end

sim_names = ["box", "zoom_inj", "zoom_dpp"]

snap_bases = ["/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000",
            "/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20e/snapdir_074/snap_074",
            "/e/ocean2/users/lboess/LocalUniverseZooms/L5/cr6p20eDpp/snapdir_074/snap_074"]

use_keys = [true, false, false]


for i in 1:length(snap_bases)
    out_dir = joinpath("/e/ocean2/users/lboess/PaperRepos/2310.13734/data/spectra", sim_names[i])
    run(snap_bases[i], use_keys[i], out_dir)
end