using SPHtoGrid
using Printf

"""
    get_slope_image(im_144, im_1440, vmin)

Computes the synchrotron slope per pixel between 144 MHz and 944 MHz.
`vmin` is an arbitrary cutoff to avoid computing pixels with spectra that are too steep to be relevant.
"""
function get_slope_image(im_144, im_1440, vmin=0.0)

    # safety cutoff
    im_144[im_144.<0.0] .= 0.0
    im_1440[im_1440.<0.0] .= 0.0

    im_144_cut = copy(im_144)
    im_144_cut[im_1440.<vmin] .= 0.0
    im_1440_cut = copy(im_1440)
    im_1440_cut[im_1440.<vmin] .= 0.0

    im_slope = (log10.(im_144_cut) .- log10.(im_1440_cut)) ./ (log10(144.e6) - log10(1.4e9))

    im_slope
end

"""
    write_slope_image(snap, orientation)
"""
function write_slope_image(map_base, snap, Bname)

    filenames = ["synch_Inu_144MHz", "synch_Inu_1.4GHz"] .* "_$Bname"

    files = [map_base * "coma_20Mpc_$snap.$filename.xz.fits"
             for filename ∈ filenames]

    println("loading data")
    im_144, par, snap_num, units = read_fits_image(files[1])
    im_1440, par, snap_num, units = read_fits_image(files[2])

    println("constructing slope image")
    im_slope = get_slope_image(im_144, im_1440, 1.e-26)

    filename = map_base * "coma_20Mpc_$snap.synch_slope_$Bname.xz.fits"
    println("saving slope image")
    write_fits_image(filename, im_slope,
        par,
        units="[]",
        snap=snap_num)

    println("done")

end

Bnames = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]
folders = "coma/" * ["box", "zoom_inj", "zoom_dpp_1e-17", "zoom_dpp_5e-17", "zoom_HB07"]  .* "/"
snaps = ["036", "012", "074", "012", "012"]


for i ∈ 5:length(folders), Bname ∈ Bnames 
    println(folders[i], " ", Bname)
    write_slope_image(map_path * folders[i], snaps[i], Bname)
end
