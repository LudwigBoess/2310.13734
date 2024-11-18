"""
    Explanatory figure for emissivity scaling with slope and p_inj
"""

include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using SpectralCRsUtility
using GadgetUnits
using Base.Threads
using PyPlot, PyPlotUtility
@info "done"

const global slope_soft = 1.e-6

"""
    find_init_norm(pressure::T, slope::T, bound_low::T, bound_up::T) where T

Find norm of first bin for a given total pressure.
"""
function find_init_norm(ϵ::T, slope::T, bound_low::T, bound_up::T) where {T}

    norm_cnst = ϵ / (4π * SpectralCRsUtility.cL * bound_low^4)

    init_norm = norm_cnst * (4 - slope) / ((bound_up / bound_low)^(4 - slope) - 1)

    if (4 - slope_soft) < slope < (4 + slope_soft)

        slope_var = (slope - 4) / slope_soft
        init_norm2 = norm_cnst / log(bound_up / bound_low)

        if !iszero(slope_var)
            return init_norm * slope_var + init_norm2 * (1 - slope_var)
        else
            return init_norm2
        end
    end

    return init_norm
end


function get_phase_map(p_range, q_range, Nbins)

    GU = GadgetPhysical()

    B = 1.0e-6
    pmax = 1.e6
    ϵ_cr = 1.0

    ϵ_cr_GU = ϵ_cr * GU.x_cgs^3 / GU.E_cgs

    j_nu = Matrix{Float64}(undef, Nbins, Nbins)

    for i = 1:Nbins, j = 1:Nbins

        p_min = p_range[i]
        q = q_range[j]

        par = CRMomentumDistributionConfig(p_min, pmax, 96)
        bounds = momentum_bin_boundaries(par)

        f_0 = find_init_norm(ϵ_cr_GU, q, p_min, pmax)
        fp = [f_0 * 10.0^(j * par.bin_width * (-q)) for j = 0:par.Nbins-1] .* GU.CR_norm

        fq = q .* ones(par.Nbins)
        cut = pmax

        j_nu[i, j] = synchrotron_emission(fp, fq, cut, B, bounds, integrate_pitch_angle=false)


    end

    return j_nu
end

function plot_scaling(q_range, p_range, j_nu, plot_name)

    c_lim = [1.e-34, 1.e-27]
    levels = 10.0 .^ LinRange(-34, -27, floor((34 - 27) * 2))

    fig = get_figure(1.0)
    plot_styling!()
    ax = gca()
    ax.set_yscale("log")
    xlabel("Spectral Slope  " * L"q")
    ylabel("Injection Momentum  " * L"\hat{p}_\mathrm{inj}" * "  [" * L"(m_e c)^{-1}" * "]")
    im = pcolormesh(q_range, p_range, j_nu,
        cmap="YlGnBu", norm=matplotlib.colors.LogNorm(vmin=c_lim[1], vmax=c_lim[2]))

    CS = contour(q_range, p_range, j_nu, levels=levels, colors="gray",
        linewidths=1, linestyles="dashed")

    get_colorbar_right(ax, im, "Synchrotron Emissivity  " * L"j_{\nu = 144 \: \mathrm{MHz}}" * " [erg s" * L"^{-1}" * " Hz" * L"^{-1}" * "cm" * L"^{-3}" * "]")

    savefig(plot_name, bbox_inches="tight", dpi=500)
    close(fig)
end

Nbins = 100
p_range = 10.0 .^ LinRange(-1, 1, Nbins)
q_range = LinRange(4, 6, Nbins)

j_nu = get_phase_map(p_range, q_range, Nbins)

plot_name = plot_path * "Fig12.pdf"

plot_scaling(q_range, p_range, j_nu, plot_name)
