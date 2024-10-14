using SPHtoGrid
using GadgetIO, GadgetUnits
using Unitful, UnitfulAstro

const global sim_path = "/gpfs/work/pn36ze/di93son/LocalUniverse/Coma/L4/cr6p16e_tricco/"
const global map_path = "/gpfs/work/pn36ze/di93son/PaperRepos/2310.13734/maps/"
fi = sim_path * "snapdir_011/snap_011"
h = read_header(fi)
c = cosmology(h)

z_coma = 0.0231
beam = 60.0u"arcsecond"

θ_beam = [beam, beam]
smooth_size = zeros(2)
for i = 1:2
    smooth_size[i] = 0.5 * ustrip(arcmin_to_kpc(c, θ_beam[i], z_coma))
end

println(smooth_size)


filenames = ["Bsim", "beta50", "01Pturb", "BFF", 
            "dyn_l", "dyn_h"]

folders = ["box"#, "zoom_inj", "zoom_dpp", "zoom_HB"
            ] .* "/"
snaps = ["000", "074", "074", "074"]

files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_Inu_144MHz_$filename.xz.fits"
         for i ∈ 1:length(folders), filename ∈ filenames]
out_files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).synch_F_beam_1'_144MHz_$filename.xz.fits"
             for i ∈ 1:length(folders), filename ∈ filenames]

for i = 1:length(files)
    println(i)
    image, par, snap, units = read_fits_image(files[i])
    map_P = synchrotron_SB_to_luminosity(image, par)
    map_F = convert_Pnu_map_to_mJy_beam(map_P, par.pixelSideLength, beam, c, z_coma) .* 100.0
    println(maximum(map_F))
    write_fits_image(out_files[i], map_F, par, units="mJy/beam"; snap)
end


"""
    CReE 
"""

using SPHtoGrid
using GadgetUnits
using Unitful, UnitfulAstro

folders = ["box"]#, "zoom_inj"]
snaps = ["000", "074"]

files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).CReE_gt1GeV.xz.fits"
         for i ∈ 1:length(folders)]

out_files = [map_path * "$(folders[i])/coma_20Mpc_$(snaps[i]).CReE_gt1GeV_L.xz.fits"
             for i ∈ 1:length(folders)]


for i = 1:length(files)
    image, par, snap, units = read_fits_image(files[i])

    image_E = surface_brightness_to_luminosity(image, par)

    write_fits_image(out_files[i], image_E, par, units="erg"; snap)
end