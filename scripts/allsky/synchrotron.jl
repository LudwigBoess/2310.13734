"""
    get_synchrotron(data, nu, Bfield_function, show_progress=false)

Calculate synchrotron emissivity for a given data set at observational frequency `nu` for a magnetic field defined by `Bfield_function`.
"""
function get_synchrotron(data, nu, Bfield_function, show_progress=false)

    Npart = length(data["CReC"])
    j_ν = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)

    if show_progress
        P = Progress(Npart)
    end

    @threads for i ∈ eachindex(j_ν)

        d = 0.0
        @inbounds for dim = 1:3
            d += (data["POS"][dim, i] - center_comov[dim])^2
        end

        d = sqrt(d) * GU.x_physical

        if d < radius_limits[1] || d > radius_limits[2]
            j_ν[i] = 0.0
            continue
        end

        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]

        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])

        B = Bfield_function(data, i)

        j_ν[i] = synchrotron_emission(norm, slope, cut, B, par, ν0=nu,
            reduce_spectrum=true,
            integrate_pitch_angle=true,
            convert_to_mJy=false)

        if show_progress
            next!(P)
        end
    end

    j_ν
end

"""
    synch_maps_of_subfile(subfile, Bfield_function, Btype)

Calculate synchrotron healpix maps for a given subfile for a magnetic field defined by `Bfield_function`.
"""
function synch_maps_of_subfile(subfile, Bfield_function, Btype)

    println("syn (B$Btype): subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    hsml = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho = read_block(snap_base * ".$subfile", "RHO", parttype=0) .* GU.rho_physical
    m = read_block(snap_base * ".$subfile", "MASS", parttype=0) .* GU.m_physical
    pos = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical

    weights = part_weight_physical(length(hsml), GU.x_cgs)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["POS", "MASS", "HSML", "RHO", "U", "BFLD", "VRMS", "CReN", "CReS", "CReC"])

    # calculate synchrotron emissivity
    j_nu = get_synchrotron(data, 144.e6, Bfield_function)
    j_nu = set_rest_to_zero(pos, j_nu)
    println("\tsynch done!\n\tmaximum = $(maximum(j_nu)) erg/s/Hz/cm^3\n\tsum = $(sum(j_nu)) erg/s/Hz/cm^3")

    map = healpix_map(pos, hsml, m, rho, j_nu, weights, show_progress=false;
        center, kernel, Nside)

    println("\tmap done!\n\tmaximum = $(maximum(map[1])) erg/s/Hz/cm^2\n\tsum = $(sum(map[1])) erg/s/Hz/cm^2")
    flush(stdout)
    flush(stderr)

    pos = hsml = rho = m = weights = nothing
    j_nu = nothing
    #P_nu = nothing
    data = nothing
    GC.gc()
    return map
end



function synch_maps_of_subfile(subfile)
    if Bfield_flag == 1
        return synch_maps_of_subfile(subfile, Bfield_sim, "sim")
    elseif Bfield_flag == 2
        return synch_maps_of_subfile(subfile, Bfield_Beta, "beta")
    elseif Bfield_flag == 3
        return synch_maps_of_subfile(subfile, Bfield_vturb, "vturb")
    elseif Bfield_flag == 4
        return synch_maps_of_subfile(subfile, Bfield_FF, "ff")
    elseif Bfield_flag == 5
        return synch_maps_of_subfile(subfile, Bfield_dyn_l, "dyn_l")
    elseif Bfield_flag == 6
        return synch_maps_of_subfile(subfile, Bfield_dyn_h, "dyn_h")
    end
end

function get_synch_filename()
    if Bfield_flag == 1
        return map_path * "allsky_synch_Pnu_144MHz_B_sim_$viewpoint.fits"
    elseif Bfield_flag == 2
        return map_path * "allsky_synch_Pnu_144MHz_B_beta50_$viewpoint.fits"
    elseif Bfield_flag == 3
        return map_path * "allsky_synch_Pnu_144MHz_B_01Pturb_$viewpoint.fits"
    elseif Bfield_flag == 4
        return map_path * "allsky_synch_Pnu_144MHz_B_FF_$viewpoint.fits"
    elseif Bfield_flag == 5
        return map_path * "allsky_synch_Pnu_144MHz_B_dyn_l_$viewpoint.fits"
    elseif Bfield_flag == 6
        return map_path * "allsky_synch_Pnu_144MHz_B_dyn_h_$viewpoint.fits"
    end
end
