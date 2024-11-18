using Base.Threads

function Bfield_sim(data, i)
    sqrt(data["BFLD"][1, i]^2 + data["BFLD"][2, i]^2 + data["BFLD"][3, i]^2)
end


function Bfield_FF(data, i)
    0.4e-6 * ∛(data["RHO"][i] * GU.rho_ncm3 * 1.e4)^2
end

function Bfield_Beta(data, i, β=50)
    γ_m1 = 5 / 3 - 1
    Pth = γ_m1 * data["RHO"][i] * data["U"][i] * GU.P_th_cgs
    return √(8π * Pth / β)
end


function Bfield_vturb(data, i)
    √(4π * data["RHO"][i] * GU.rho_cgs) * data["VRMS"][i] * GU.v_cgs
end

function B_fit(x)
    p = [-16.37664752379627, -16.001138655636378, -8.072683572197233, -1.7120945476093201, -0.1333534527683321]
    return (p[1] + p[2] * x + p[3] * x^2 + p[4] * x^3 + p[5] * x^4)
end

function Bfield_dyn_l(data, i)
    rho = data["RHO"][i] * GU.rho_ncm3
    if rho < 1.e-4        
       return 10.0^B_fit(log10(rho))
    else
        return 2.5e-6 * √(rho * 1e3)
    end
end

function Bfield_dyn_h(data, i)
    rho = data["RHO"][i] * GU.rho_ncm3
    2.5e-6 * √(rho * 1e3)
end

function run_Bfld_map_of_subfile(subfile, blocks, Bfld_function, Btype)

    println("B$Btype: subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    hsml = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho = read_block(snap_base * ".$subfile", "RHO", parttype=0) .* GU.rho_physical
    m = read_block(snap_base * ".$subfile", "MASS", parttype=0) .* GU.m_physical
    pos = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0) for block ∈ blocks)
    B = Vector{Float64}(undef, length(rho))

    @threads for i ∈ eachindex(rho)
        B[i] = Bfld_function(data, i)
    end
    B = set_rest_to_zero(pos, B)
    println("\tB done!\tmaximum = $(maximum(B)) muG")

    map = healpix_map(pos, hsml, m, rho, B, rho, show_progress=false;
        center, kernel, Nside)

    sel = findall(map[2] .> 0.0)
    println("\tmap done!\tmaximum = $(maximum(map[1][sel]./map[2][sel])) G")
    println("\tBefore GC: Available Memory: $(Sys.free_memory() / 2^20) MB -> $( Sys.free_memory() / Sys.total_memory() * 100) %")

    pos = hsml = rho = m = nothing
    B = nothing
    data = nothing
    GC.gc()
    println("\tAfter GC: Available Memory: $(Sys.free_memory() / 2^20) MB -> $( Sys.free_memory() / Sys.total_memory() * 100) %")

    return map
end


function Bfld_map_of_subfile(subfile)
    if Bfield_flag == 1
        return run_Bfld_map_of_subfile(subfile, ["BFLD"], Bfield_sim, "sim")
    elseif Bfield_flag == 2
        return run_Bfld_map_of_subfile(subfile, ["RHO", "U"], Bfield_Beta, "beta")
    elseif Bfield_flag == 3
        return run_Bfld_map_of_subfile(subfile, ["RHO", "VRMS"], Bfield_vturb, "vturb")
    elseif Bfield_flag == 4
        return run_Bfld_map_of_subfile(subfile, ["RHO"], Bfield_FF, "ff")
    elseif Bfield_flag == 5
        return run_Bfld_map_of_subfile(subfile, ["RHO"], Bfield_dyn_l, "dyn_l")
    elseif Bfield_flag == 6
        return run_Bfld_map_of_subfile(subfile, ["RHO"], Bfield_dyn_h, "dyn_h")
    end
end


function get_B_filename()
    if Bfield_flag == 1
        return map_path * "allsky_B_sim_$viewpoint.fits"
    elseif Bfield_flag == 2
        return map_path * "allsky_B_beta50_$viewpoint.fits"
    elseif Bfield_flag == 3
        return map_path * "allsky_B_Pturb_$viewpoint.fits"
    elseif Bfield_flag == 4
        return map_path * "allsky_B_FF_$viewpoint.fits"
    elseif Bfield_flag == 5
        return map_path * "allsky_B_dyn_l_$viewpoint.fits"
    elseif Bfield_flag == 6
        return map_path * "allsky_B_dyn_h_$viewpoint.fits"
    end
end