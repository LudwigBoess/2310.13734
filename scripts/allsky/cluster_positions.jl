using GadgetIO, GadgetUnits
using Healpix
using SPHtoGrid, SPHKernels


"""
    New
"""
named_clusters_pos = [
    ClusterPositions("Centaurus", [225160.875, 256677.5625, 242031.765625], 1771.7005615234375),
    ClusterPositions("Coma", [245040.453125, 327781.84375, 246168.6875], 2177.5625),
    ClusterPositions("Fornax", [252454.0, 238890.296875, 233772.890625], 714.1318359375),
    ClusterPositions("Norma", [202408.359375, 245067.234375, 249614.09375], 798.90478515625),
    ClusterPositions("Ophiuchus", [176696.22, 266631.8, 315976.38], 1067.0657),
    ClusterPositions("Perseus", [307184.78125, 247627.078125, 230736.484375], 1813.0419921875),
    ClusterPositions("Virgo", [244450.578125, 255851.78125, 253668.1875], 1711.311767578125), 
    ClusterPositions("A119", [321548.9375, 128183.8125, 230358.765625], 1775.019287109375),
    ClusterPositions("A539", [302446.71875, 254504.890625, 173170.796875], 1522.255859375),
    ClusterPositions("A569", [298514.34, 278775.4, 226338.19], 800.08154),
    ClusterPositions("A576", [379893.44, 313863.3, 190307.94], 1586.8203),
    ClusterPositions("A634", [289826.4, 290567.34, 222990.25], 879.5112),
    ClusterPositions("A754", [217824.12, 338906.38, 113371.805], 846.7364),
    ClusterPositions("A779", [291807.2, 304191.44, 207839.58], 675.18677),
    ClusterPositions("A1060", [184219.84, 260534.86, 204199.1], 1452.2032),
    ClusterPositions("A1185", [272634.6, 357607.56, 228079.97], 1351.1272),
    ClusterPositions("A1228", [278397.38, 370507.72, 234770.06], 1169.6703),
    ClusterPositions("A1257", [278496.03, 366316.66, 241958.06], 957.47766),
    ClusterPositions("A1267", [248213.36, 357303.75, 234187.6], 948.92285),
    ClusterPositions("A1367", [245088.19, 316851.6, 236572.25], 1086.3784),
    ClusterPositions("A1631a", [135390.72, 379839.88, 251627.33], 1040.5299),
    ClusterPositions("A1644", [141584.67, 370061.25, 253373.22], 1503.0912),
    ClusterPositions("A1736", [121316.18, 340799.75, 253865.84], 1208.1821),
    ClusterPositions("A2040", [187733.44, 345981.53, 343366.4], 1214.6323),
    ClusterPositions("A2063", [191684.47, 309146.44, 329393.6], 1041.5564),
    ClusterPositions("A2256", [363955.6, 341751.94, 334906.94], 1233.0735),
    ClusterPositions("A2107", [202172.78, 333935.25, 355058.25], 1013.9514),
    ClusterPositions("A2147", [219706.1, 316873.25, 338450.2], 1109.9745),
    ClusterPositions("A2162", [248152.6, 300697.8, 330477.4], 1003.3834),
    ClusterPositions("A2197", [262659.56, 335529.16, 313575.25], 931.2115),
    ClusterPositions("A2199", [264176.1, 325519.38, 312852.2], 1062.0698),
    ClusterPositions("A2319", [326681.03, 277891.8, 392259.53], 1057.7314),
    ClusterPositions("A2877", [239844.34375, 169326.65625, 254522.625], 1082.710693359375),
    ClusterPositions("A3532", [102092.87, 364356.22, 236999.48], 1001.11035),
    ClusterPositions("A3566", [118897.08, 341342.12, 259078.03], 927.4432),
    ClusterPositions("A3570", [150716.69, 310394.5, 249586.64], 1744.2521),
    ClusterPositions("A3571", [141463.22, 302101.4, 256304.31], 1126.0339),
    ClusterPositions("A3581", [164433.22, 276585.56, 243083.98], 864.4208), 
    ClusterPositions("RXJ1030.0-3521", [212907.31, 259784.52, 232045.38], 1171.2385),
    ClusterPositions("RXJ1253.0-0912", [208313.88, 291612.12, 257487.16], 983.10406),
    ClusterPositions("RXJ1315.3-1623", [210838.38, 263100.94, 245662.97], 706.7438),
    ClusterPositions("RXJ1334.3+3441", [262976.2, 324173.12, 262512.78], 878.4004),
    ClusterPositions("RXJ1740.5+3539", [249403.25, 310609.78, 390071.34], 1804.4414), 
    ClusterPositions("AWM_5", [234972.8, 294223.66, 352418.16], 1389.4089),
    ClusterPositions("MKW_8", [218508.69, 302728.25, 294728.9], 1025.9329), 
    ClusterPositions("Z4803", [220816.84, 297480.06, 238753.45], 810.2528)
]

function get_cluster_outlines(clusters, GU)


    pos  = Matrix{Float64}(undef, 3, length(clusters))
    hsml = Vector{Float64}(undef, length(clusters))

    for i = 1:length(clusters) 
        pos[:,i] = clusters[i].pos  .* GU.x_physical
        hsml[i]  = clusters[i].rvir * GU.x_physical
    end

    weights  = ones(length(hsml))
    center = [246.980, 245.478, 254.286] .* 1.e3 .* GU.x_physical

    kernel = Tophat()
    radius_limits = [1_000.0, Inf]

    allsky_map, weight_map = healpix_map(pos, hsml, weights, weights, show_progress=true; 
                            center, kernel, Nside, radius_limits)

    @inbounds for i ∈ eachindex(allsky_map)
        if !isnan(weight_map[i]) && !iszero(weight_map[i]) && !isinf(weight_map[i])
            allsky_map[i]  /= weight_map[i]
        end
    end

    return allsky_map
end


map = get_cluster_outlines(named_clusters_pos, GU)
