include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "shared.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Statistics
using Printf
using ProgressMeter
using PyCall
# needs to by imported by hand to make inset axis
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
axes_divider = pyimport("mpl_toolkits.axes_grid1.axes_divider")
@info "done"

function Bfield_FF(rho)
    0.1e-6 * ∛(rho * 1.e4)^2 * 4
end


function Bfield_dyn_h(rho)
    2.5e-6 * √(rho * 1e3)
end

function B_fit(x)
    p = [-16.37664752379627, -16.001138655636378, -8.072683572197233, -1.7120945476093201, -0.1333534527683321]
    return (p[1] + p[2] * x + p[3] * x^2 + p[4] * x^3 + p[5] * x^4)
end

function Bfield_dyn_l(rho)
    if rho < 1.e-4        
       return 10.0^B_fit(log10(rho))
    else
        return 2.5e-6 * √(rho * 1e3)
    end
end


function plot_models(phase_map_path, plot_name)

    Bfield_models = [L"B_\mathrm{sim}",
        L"B_{\beta}",
        L"B_{\mathcal{F}}",
        L"B_\mathrm{ff}",
        L"B_\mathrm{dyn ↓}",
        L"B_\mathrm{dyn ↑}"]

    Bfield_functions = [Bfield_FF, Bfield_dyn_l, Bfield_dyn_h]

    filename = phase_map_path .* ["box/bin_1D_$B.dat"
                             for B ∈ ["B_sim", "B_beta50", "B_01Pturb"] ]

    sm2 = plt.cm.ScalarMappable(cmap=PyPlot.cm.magma,
        norm=plt.Normalize(vmin=0, vmax=6.5))
    sm2.set_array([])

    lw = 3
    alpha_ref = 1.0

    fig = get_figure(1.0)
    plot_styling!()
        ax = gca()

        axis_ticks_styling!(ax)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim([1.e-7, 1.0e-2])
        ax.set_ylim([1.e-14, 1.e-4])
        ax.xaxis.set_major_locator(plt.LogLocator(10))

        ax.set_facecolor("white")
        ax.set_xlabel("Electron Density  " * L"n_e" * "  [cm" * L"^{-3}]")
        ax.set_ylabel("Magnetic Field Strength  " * L"B" * "  [G]")

        for i = 1:3
            ne_bins, B_mean = read_1D_data(filename[i])
            sel = findall(ne_bins .< 5.e-9)
            B_mean[sel] .= NaN
            sel = findall(ne_bins .> 3.e-1)
            B_mean[sel] .= NaN

            ax.plot(ne_bins, B_mean, alpha=alpha_ref,
                color=sm2.to_rgba(i), lw=lw, label=Bfield_models[i])
        end

        for i = 1:3
            ne_bins, B_mean = read_1D_data(filename[i])
            B = @. Bfield_functions[i](ne_bins)
            ax.plot(ne_bins, B, alpha=alpha_ref,
                color=sm2.to_rgba(i+3), lw=lw, label=Bfield_models[i+3])
        end

    legend(frameon=false, loc="lower right")

    savefig(plot_name, bbox_inches="tight", transparent=false)
    close(fig)

end

phase_map_path = data_path * "phase_maps/"
plot_name = plot_path * "FigD1.pdf"

plot_models(phase_map_path, plot_name)