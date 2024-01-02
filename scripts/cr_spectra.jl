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
    spec::CRMomentumDistribution
    j_nu::Vector{Float64}
end

function read_spectra()

    filename = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/spectra.dat"
    f = open(filename, "r")
    Nids = read(f, Int64)
    ids = read!(f, Vector{UInt64}(undef, Nids))
    Nfiles = read(f, Int64)
    spectra = Matrix{CR}(undef, Nids, Nfiles)

    for Nfile ∈ 1:Nfiles, Nid ∈ 1:Nids
        bound = read!(f, Vector{Float64}(undef, 49))
        norm = read!(f, Vector{Float64}(undef, 48))
        j_nu = read!(f, Vector{Float64}(undef, 50))

        spectra[Nid, Nfile] = CR(CRMomentumDistribution(bound, norm), j_nu)
    end
    close(f)

    return ids, spectra
end

function plot_spectra(spectra, ids, t)

    sm = plt.cm.ScalarMappable(cmap=PyPlot.cm.jet,
        norm=plt.Normalize(vmin=7.0, vmax=14.0))
    sm.set_array([])

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

        @showprogress for i ∈ 1:length(t)
            ax.plot(spectra[Nid, i].spec.bound[1:end-1], spectra[Nid, i].spec.norm, 
                    c=sm.to_rgba(t[i]))
        end

        #plot([1.e4, 1.e5], [1.e-29, 1.e-29 * 10.0^(-q0)], color="k", linestyle="--", linewidth=2)
        #text(2.e4, 1.e-31, L"q_0 \: = " * "$(-q0)", rotation=-45)


        subplot(get_gs(gs, 0, 1))
        ax = gca()
        ax.set_xlim([1.e7, 1.e10])
        ax.set_ylim([1.e-51, 1.e-36])
        ax.set_xscale("log")
        ax.set_yscale("log")
        axis_ticks_styling!(ax)

        xlabel("Obs. Frequency  " * L"\nu" * " [Hz]")
        ylabel("Synch. Emissivity  " * L"j_\nu" * " [erg s" * L"^{-1}" * "Hz" * L"^{-1}" * "cm" * L"^{-3}" * "]")
        ν_arr = 10.0 .^ (LinRange(7, 10, 50))
        @showprogress for i ∈ 1:length(t)
            ax.plot(ν_arr, spectra[Nid, i].j_nu, c=sm.to_rgba(t[i]))
        end

        #plot([1.e9, 1.e10], [2.e-27, 2.e-27 * 10.0^(-α(q0))], color="k", linestyle="--", linewidth=2)
        #text(2.e9, 1.e-27, L"\alpha_0 \: = " * "$(-α(q0))", rotation=-38)
        subplot(get_gs(gs, 0, 2))
        cax = gca()
        cb = colorbar(sm, cax=cax, fraction=0.046)
        cb.set_label("Time  " * L"t" * " [Gyr]")
        cb.ax.tick_params(
            direction="in",
            which="major",
            size=6, width=1
        )

        plot_name = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/Plots/spectra/spec_$(ids[Nid]).png"
        savefig(plot_name, bbox_inches="tight")
        close(fig)
    end

end

ids, spectra = read_spectra()
times = readdlm("/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/data/times.txt")[:,1]

plot_spectra(spectra, ids, times)