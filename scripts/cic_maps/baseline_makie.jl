using CairoMakie
using Printf
using SPHtoGrid
using ColorSchemes
using LaTeXStrings

function plot_relic_evolution(plot_name)

    @info "plotting..."

    Bfield_names = ["Bsim", "beta50", "01Pturb", "BFF", "dyn_l", "dyn_h"]
    files = @. "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/coma/zoom_dpp_5e-17/coma_20Mpc_$(@sprintf("%03i", 12)).baseline_" * Bfield_names * ".xz.fits"

    Ncols = 6
    Nrows = 1

    f = Figure()#size = (800, 1200))
    cb = GridLayout(f[1, 1], 1, Ncols, alignmode=Inside())
    ax = GridLayout(f[2, 1], Nrows, Ncols, alignmode=Inside())

    _hm = nothing # dummy variable
    for row = 1:Nrows, col = 1:1
        map, par, snap, dummy = read_fits_image(files[row])
        sel = findall(.!isinf.(map) .&& .!isnan.(map))
        println(maximum(map[sel]))
        println(minimum(map[sel]))


        ca = Axis(ax[row, col], aspect = DataAspect(), 
            bottomspinecolor=:white, topspinecolor=:white, leftspinecolor=:white, rightspinecolor=:white)
        hidedecorations!(ca, grid=false, label=false)
        _hm = heatmap!(ca, map, colormap=ColorSchemes.diverging_bky_60_10_c30_n256, 
                colorrange=(-1.0e-17, 1.0e-17),
                colorscale=Makie.Symlog10(1.e-19),
            margin=(0.0, 0.0))
    end


    for col = 1:Ncols
        colsize!(ax, col, Aspect(1, 1.0))
    end

    println(_hm)

    colgap!(ax, 0)
    rowgap!(ax, 0)
    Colorbar(cb[1,:], _hm, vertical=false,
        minorticksvisible=true, minortickalign=1, tickalign=1,
        label="Rel. div B")#, font="Times New Roman")


    #ref_axis = ax.scene.px_area[].widths[2]
    #colsize!(f.layout, 1, ref_axis.scene.px_area[].widths[2])

    colsize!(cb, 1, Aspect(1, 1))
    # for row = 2:Nrows+1
    #     rowsize!(ax, row, Aspect(1, 1.0))
    # end

    #colsize!(Axis(f[1, 1]), 1, Aspect(1, 1.0))


    resize_to_layout!(f)

    save(plot_name, f)
end


plot_name = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/baseline.png"

plot_relic_evolution(plot_name)

