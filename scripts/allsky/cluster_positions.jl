using GadgetIO, GadgetUnits
using Healpix
using SPHtoGrid, SPHKernels
using DelimitedFiles
using Printf
using ProgressMeter
using SkyCoords

# simulation settings
const global sim_path = "/gpfs/work/pn68va/di67meg/LocalUniverse/"
const global snap = 36
const global snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"

# general map settings
const global Nside = 2048
const global Nfiles = 2048
const global map_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/allsky/"
const global viewpoint = "slow_1"

const global h = GadgetIO.read_header(snap_base)
const global GU = GadgetPhysical(h)
const global kpc2cm = 3.085678e21


# SLOW 1 paper
const global center_comov = [247.980, 245.480, 255.290] .* 1.e3
const global center = center_comov .* GU.x_physical

function galactic_pixel(res, i)

    b, l = pix2angRing(res, i)
    supergal = SuperGalCoords(l, colat2lat(b))
    gal = convert(GalCoords, supergal)

    return ang2pixRing(res, π / 2 - gal.b, gal.l)

end

function supergalactic_angle(res, i)

    b, l = pix2angRing(res, i)
    gal = GalCoords(l, colat2lat(b))
    supergal = convert(SuperGalCoords, gal)

    return (π / 2 - supergal.b), supergal.l
end


function rotate_healpix_map(m)

    Nside = npix2nside(length(m[:]))

    res = Healpix.Resolution(Nside)

    m_supergal = HealpixMap{Float64,RingOrder}(Nside)

    @showprogress for i ∈ 1:nside2npix(Nside)
        θ, ϕ = supergalactic_angle(res, i)
        m_supergal[i] = interpolate(m, θ, ϕ)
    end

    return m_supergal
end


function make_healpix_map(pos, hsml, center, Nside)

    kernel = Tophat()
    weights = ones(length(hsml))

    rho = @. weights / (4π / 3 * hsml^3)

    bin_q = ones(length(hsml))
    bin_q .= Inf

    Pos = copy(pos)

    @info "making map"

    allsky_map, weight_map = healpix_map(Pos, hsml, weights, rho, bin_q, weights, show_progress=true;
        center, kernel, Nside)

    @inbounds for i ∈ eachindex(allsky_map)
        if !isnan(weight_map[i]) && !iszero(weight_map[i]) && !isinf(weight_map[i])
            allsky_map[i] /= weight_map[i]
        end
    end

    @inbounds for i ∈ eachindex(allsky_map)
        if isinf(allsky_map[i])
            allsky_map[i] = 1.0
        end
    end

    allsky_filename = "/gpfs/work/pn68va/di67meg/Paper/LocalUniverse/synchrotron_web/maps/allsky_slow1/cluster_contours_r500_slow_1.fits"
    if isfile(allsky_filename)
        rm(allsky_filename)
    end
    saveToFITS(allsky_map, allsky_filename)

    @info "rotating"
    m = rotate_healpix_map(allsky_map)
    #m = allsky_map

    allsky_filename = "/gpfs/work/pn68va/di67meg/Paper/LocalUniverse/synchrotron_web/maps/allsky_slow1/cluster_contours_r500_slow_1_gal.fits"
    if isfile(allsky_filename)
        rm(allsky_filename)
    end
    saveToFITS(m, allsky_filename)

    @info "projecting"
    cluster_proj, mask, maskflag = project(mollweideprojinv, m, 2Nside, Nside)
    return cluster_proj
end


fi = "/gpfs/work/pn68va/di67meg/Paper/LocalUniverse/synchrotron_web/data/data_elly.txt"

data = readdlm(fi, skipstart=1)
pos = copy(transpose(data[:, 3:5])) .* GU.x_physical .|> Float64
hsml = data[:, 6] .* GU.x_physical .|> Float64
contour_image = make_healpix_map(pos, hsml, center, 2048)
