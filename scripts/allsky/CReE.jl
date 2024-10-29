function get_CReE_seed(filename)

    data = Dict(block => read_block(filename, block, parttype=0)
                for block ∈ ["MASS", "RHO", "CReN", "CReS", "CReC"])

    h = read_header(filename)
    GU = GadgetPhysical(h)

    Npart = length(data["CReC"])
    CReE = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(1.0, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    @threads for i = 1:Npart

        norm = 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])
        rho = data["RHO"][i] * GU.rho_physical

        CReE[i] = cr_energy_in_range(norm, slope, cut, rho, bounds,
            Emin=5.e-3, Emax=0.5)

    end

    # convert to erg
    return CReE .* data["RHO"] .* GU.E_cgs ./ GU.x_cgs^3
end

function CReE_seed_maps_of_subfile(subfile)

    println("CReP_seed: subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    hsml = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho = read_block(snap_base * ".$subfile", "RHO", parttype=0) .* GU.rho_physical
    m = read_block(snap_base * ".$subfile", "MASS", parttype=0) .* GU.m_physical

    pos = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical

    CReP = get_CReE_seed(snap_base * ".$subfile")
    CReP = set_rest_to_zero(pos, CReP)

    map = healpix_map(pos, hsml, m, rho, CReP, rho, show_progress=false, calc_mean=false;
        center, kernel, Nside)

    pos = hsml = rho = m = nothing
    CReP = nothing

    GC.gc()
    return map
end


function get_CReE(filename)

    data = Dict(block => read_block(filename, block, parttype=0)
                for block ∈ ["MASS", "RHO", "CReN", "CReS", "CReC"])

    h = read_header(filename)
    GU = GadgetPhysical(h)

    Npart = length(data["CReC"])
    CReE = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(1.0, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    @threads for i = 1:Npart

        norm = 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])
        #rho    = data["RHO"][i] * GU.rho_physical

        CReE[i] = cr_energy_in_range(norm, slope, cut, 1.0, bounds)

    end

    # convert to erg
    return CReE .* (GU.E_cgs / GU.x_cgs^3)

end

function CReE_maps_of_subfile(subfile)

    println("CReE: subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    hsml = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho = read_block(snap_base * ".$subfile", "RHO", parttype=0) .* GU.rho_physical
    m = read_block(snap_base * ".$subfile", "MASS", parttype=0) .* GU.m_physical
    pos = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical

    CReE = get_CReE(snap_base * ".$subfile")
    CReE = set_rest_to_zero(pos, CReE)

    weights = part_weight_physical(length(hsml), GU.x_cgs)

    map = healpix_map(pos, hsml, m, rho, CReE, weights, show_progress=false, calc_mean=false;
        center, kernel, Nside)

    # force garbage collect
    pos = hsml = rho = m = nothing
    CReE = nothing
    GC.gc()

    return map
end
