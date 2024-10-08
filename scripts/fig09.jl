include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"

using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using SpectralCRsUtility
using Printf
using PyCall
using DelimitedFiles
using ProgressMeter
cm = pyimport("cmasher")
@info "done"

struct CR
    rho::Float64
    B::Float64
    Mach::Float64
    spec::CRMomentumDistribution
    j_nu::Vector{Float64}
end

function read_spectra()

    filename = data_path * "spectra_073.dat"
    f = open(filename, "r")
    Nids = read(f, Int64)
    ids = read!(f, Vector{UInt64}(undef, Nids))
    Nfiles = read(f, Int64)
    spectra = Matrix{CR}(undef, Nids, Nfiles)

    for Nfile ∈ 1:Nfiles, Nid ∈ 1:Nids
        rho = read(f, Float64)
        B = read(f, Float64)
        Mach = read(f, Float64)

        bound = read!(f, Vector{Float64}(undef, 41))
        norm = read!(f, Vector{Float64}(undef, 40))
        j_nu = read!(f, Vector{Float64}(undef, 50))

        spectra[Nid, Nfile] = CR(rho, B, Mach, CRMomentumDistribution(bound, norm), j_nu)
    end
    close(f)

    return ids, spectra
end

function plot_spectra_2x1(spectra, ids, t)

    sm = plt.cm.ScalarMappable(cmap=cm.sapphire,
        norm=plt.Normalize(vmin=7.0, vmax=14.0))
    sm.set_array([])

    bounds = momentum_bin_boundaries(CRMomentumDistributionConfig(0.1, 1.e5, 24))

    for Nid ∈ 1:length(ids)

        println("id $(ids[Nid])")
        fig = get_figure(2.5)
        plot_styling!()
        gs = plt.GridSpec(1, 3, hspace=0.0, wspace=0.3, width_ratios=[1, 1, 0.05], figure=fig)

        subplot(get_gs(gs, 0, 0))
        ax = gca()
        ax.set_xlim([0.8e2, 1.2e5])
        ax.set_ylim([1.e-33, 1.e-21])
        ax.set_xscale("log")
        ax.set_yscale("log")
        axis_ticks_styling!(ax)

        xlabel("Dimensionless Momentum  " * L"\hat{p}" * " [" * L"(m_e \: c)^{-1}" * "]")
        ylabel("Distribution Function  " * L"f(\hat{p})" * " [arb. units]")

        for bound ∈ bounds 
            axvline(bound, color="k", linestyle="--", linewidth=1, alpha=0.2)
        end

        @showprogress for i ∈ 1:length(t) #2:2:length(t)-2
            ax.plot(spectra[Nid, i].spec.bound[1:end-1], spectra[Nid, i].spec.norm, 
                    c=sm.to_rgba(t[i]))
        end

        plot([1.e3, 1.e4], [1.e-25, 1.e-25 * 10.0^(-4)], color="k", linestyle="--", linewidth=2)
        text(3.e3, 1.e-27, L"q_0 \: = " * "$(-4)", rotation=-45)

        get_cr_energy_axis!(ax, "e")

        subplot(get_gs(gs, 0, 1))
        ax = gca()
        ax.set_xlim([1.e7, 1.e10])
#        ax.set_ylim([1.e-41, 1.e-37])
        #ax.set_ylim([1.e-43, 1.e-39])

        ax.set_xscale("log")
        ax.set_yscale("log")
        axis_ticks_styling!(ax)

        xlabel("Obs. Frequency  " * L"\nu" * " [Hz]")
        ylabel("Synch. Emissivity  " * L"j_\nu" * " [erg s" * L"^{-1}" * "Hz" * L"^{-1}" * "cm" * L"^{-3}" * "]")
        
        ν_arr = 10.0 .^ (LinRange(7, 10, 50))
        @showprogress for i ∈ 2:2:length(t)-2
            ax.plot(ν_arr, spectra[Nid, i].j_nu, c=sm.to_rgba(t[i]))
        end

        plot([1.e8, 1.e9], [2.e-40, 2.e-40 * 10.0^(-0.5)], color="k", linestyle="--", linewidth=2)
        text(2.e8, 1.e-40, L"\alpha_0 \: = " * "$(-0.5)", rotation=-20)


        subplot(get_gs(gs, 0, 2))
        cax = gca()
        cb = colorbar(sm, cax=cax, fraction=0.046)
        cb.set_label("Time  " * L"t" * " [Gyr]")
        cb.ax.tick_params(
            direction="in",
            which="major",
            size=6, width=1
        )

        #plot_name = plot_path * "Fig09.pdf"
        plot_name = plot_path * "Fig09_$(ids[Nid]).png"
        savefig(plot_name, bbox_inches="tight")
        close(fig)
    end

end


ids, spectra = read_spectra()
times = readdlm(data_path * "times.txt")[:,1]

plot_spectra_2x1(spectra, ids, times)