
"""
    bin_2D( x_q, y_q, x_lim, y_lim, bin_q=nothing; 
            calc_mean::Bool=true, Nbins::Int=100,
            show_progress::Bool=true)

Get a 2D histogram of `x_q` and `y_q` in the limits `x_lim`, `y_lim` over a number if bins `Nbins`.
If you want to bin an additional quantity `bin_q` you can pass it as an optional argument. 
If `calc_mean` is set to true it will calculate the mean of the quantity per bin, otherwise the sum is returned.
"""
function bin_2D(x_q, y_q, x_lim, y_lim, bin_q=nothing;
    calc_mean::Bool=true, Nbins::Int=100,
    show_progress::Bool=true)

    # get bin spacing
    dx = (x_lim[2] - x_lim[1]) / Nbins
    dy = (y_lim[2] - y_lim[1]) / Nbins

    # allocate Nbins x Nbins matrix filled with zeros
    phase_map_count = zeros(Int64, Nbins, Nbins)

    if !isnothing(bin_q)
        # if a quantity should be mapped allocate Nbins x Nbins matrix filled with zeros
        phase_map = zeros(eltype(x_q[1]), Nbins, Nbins)
    end

    # optional progress meter
    if show_progress
        P = Progress(size(x_q, 1))
        idx_p = 0
    end

    @inbounds for i = 1:size(x_q, 1)

        x_bin = 1 + floor(Int64, (x_q[i] - x_lim[1]) / dx)
        y_bin = 1 + floor(Int64, (y_q[i] - y_lim[1]) / dy)

        if (1 <= x_bin <= Nbins) && (1 <= y_bin <= Nbins)

            phase_map_count[y_bin, x_bin] += 1

            if !isnothing(bin_q)
                # sum up to total binned quantity
                phase_map[y_bin, x_bin] += bin_q[i]
            end
        end

        # update progress meter
        if show_progress
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
        end
    end

    if isnothing(bin_q)
        # return simple histogram
        return phase_map_count
    else
        # if mean of binned quantity should be calculated
        if calc_mean
            # reduce map
            @inbounds for i = 1:Nbins, j = 1:Nbins
                phase_map[j, i] /= phase_map_count[j, i]
            end
        end

        # return binned quantity
        return phase_map
    end
end

"""
    bin_2D!( phase_map_count,  
             x_q, y_q, x_lim, y_lim; 
             Nbins::Int=100,
             show_progress::Bool=true)

Get a 2D histogram of `x_q` and `y_q` in the limits `x_lim`, `y_lim` over a number if bins `Nbins` for pre-allocated array `phase_map_count`.
Should be used if computing a 2D phase map over multiple files.
"""
function bin_2D!(phase_map_count,
    x_q, y_q, x_lim, y_lim;
    Nbins::Int=100,
    show_progress::Bool=true)

    # get bin spacing
    dx = (x_lim[2] - x_lim[1]) / Nbins
    dy = (y_lim[2] - y_lim[1]) / Nbins

    # optional progress meter
    if show_progress
        P = Progress(size(x_q, 1))
        idx_p = 0
    end

    @inbounds for i = 1:size(x_q, 1)

        x_bin = 1 + floor(Int64, (x_q[i] - x_lim[1]) / dx)
        y_bin = 1 + floor(Int64, (y_q[i] - y_lim[1]) / dy)

        if (1 <= x_bin <= Nbins) && (1 <= y_bin <= Nbins)

            phase_map_count[y_bin, x_bin] += 1
        end

        # update progress meter
        if show_progress
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
        end
    end

    return phase_map_count

end


"""
    bin_2D_quantity!( phase_map_count, phase_map, 
                      x_q, y_q, bin_q, x_lim, y_lim; 
                      Nbins::Int=100,
                      show_progress::Bool=true )

Get a 2D histogram of `x_q` and `y_q` in the limits `x_lim`, `y_lim` over a number of bins `Nbins` for pre-allocated arrays `phase_map_count` and `phase_map`.
Should be used if computing a 2D phase map over multiple files.
"""
function bin_2D_quantity!(phase_map_count, phase_map,
    x_q, y_q, bin_q, x_lim, y_lim;
    Nbins::Int=100, show_progress::Bool=true)

    # get logarithmic bin spacing
    dx = (x_lim[2] - x_lim[1]) / Nbins
    dy = (y_lim[2] - y_lim[1]) / Nbins

    # optional progress meter
    if show_progress
        P = Progress(size(x_q, 1))
        idx_p = 0
    end

    @inbounds for i = 1:size(x_q, 1)

        x_bin = 1 + floor(Int64, (x_q[i] - x_lim[1]) / dx)
        y_bin = 1 + floor(Int64, (y_q[i] - y_lim[1]) / dy)

        if (1 <= x_bin <= Nbins) && (1 <= y_bin <= Nbins)

            phase_map_count[y_bin, x_bin] += 1

            phase_map[y_bin, x_bin] += bin_q[i]

        end

        # update progress meter
        if show_progress
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
        end
    end

    return phase_map_count, phase_map

end





"""
    bin_2D_log( x_q, y_q, x_lim, y_lim, bin_q=nothing; 
            calc_mean::Bool=true, Nbins::Int=100,
            show_progress::Bool=true)

Get a 2D histogram of `x_q` and `y_q` in the limits `x_lim`, `y_lim` over a number if bins `Nbins` in ``log``-space.
If you want to bin an additional quantity `bin_q` you can pass it as an optional argument. 
If `calc_mean` is set to true it will calculate the mean of the quantity per bin, otherwise the sum is returned.
"""
function bin_2D_log(x_q, y_q, x_lim, y_lim, bin_q=nothing;
    calc_mean::Bool=true, Nbins::Int=100,
    show_progress::Bool=true)

    # get logarithmic bin spacing
    dlogx = (log10(x_lim[2]) - log10(x_lim[1])) / Nbins
    dlogy = (log10(y_lim[2]) - log10(y_lim[1])) / Nbins

    # allocate Nbins x Nbins matrix filled with zeros
    phase_map_count = zeros(Int64, Nbins, Nbins)

    if !isnothing(bin_q)
        # if a quantity should be mapped allocate Nbins x Nbins matrix filled with zeros
        phase_map = zeros(eltype(x_q[1]), Nbins, Nbins)
    end

    # optional progress meter
    if show_progress
        P = Progress(size(x_q, 1))
        idx_p = 0
    end

    # loop over all entries
    @inbounds for i = 1:size(x_q, 1)

        # floor division to get relevant bin
        # Julia is 1-based -> additional 1
        x_bin = 1 + floor(Int64, (log10(x_q[i]) - log10(x_lim[1])) / dlogx)
        y_bin = 1 + floor(Int64, (log10(y_q[i]) - log10(y_lim[1])) / dlogy)

        # check if bin is in range to check if particle is relevant
        if (1 <= x_bin <= Nbins) && (1 <= y_bin <= Nbins)

            # count up histogram storage
            phase_map_count[y_bin, x_bin] += 1

            if !isnothing(bin_q)
                # sum up to total binned quantity
                phase_map[y_bin, x_bin] += bin_q[i]
            end
        end

        # update progress meter if defined
        if show_progress
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
        end
    end

    if isnothing(bin_q)
        # return simple histogram
        return phase_map_count
    else
        # if mean of binned quantity should be calculated
        if calc_mean
            # reduce map
            @inbounds for i = 1:Nbins, j = 1:Nbins
                phase_map[j, i] /= phase_map_count[j, i]
            end
        end

        # return binned quantity
        return phase_map
    end
end

"""
    bin_2D_log!( phase_map_count,  
                 x_q, y_q, x_lim, y_lim; 
                 Nbins::Int=100,
                show_progress::Bool=true)

Get a 2D histogram of `x_q` and `y_q` in the limits `x_lim`, `y_lim` over a number if bins `Nbins` in ``log``-space for pre-allocated array `phase_map_count`.
Should be used if computing a 2D phase map over multiple files.
"""
function bin_2D_log!(phase_map_count,
    x_q, y_q, x_lim, y_lim;
    Nbins::Int=100,
    show_progress::Bool=true)

    # get logarithmic bin spacing
    dlogx = (log10(x_lim[2]) - log10(x_lim[1])) / Nbins
    dlogy = (log10(y_lim[2]) - log10(y_lim[1])) / Nbins

    # optional progress meter
    if show_progress
        P = Progress(size(x_q, 1))
        idx_p = 0
    end

    @inbounds for i = 1:size(x_q, 1)

        x_bin = 1 + floor(Int64, (log10(x_q[i]) - log10(x_lim[1])) / dlogx)
        y_bin = 1 + floor(Int64, (log10(y_q[i]) - log10(y_lim[1])) / dlogy)

        if (1 <= x_bin <= Nbins) && (1 <= y_bin <= Nbins)

            phase_map_count[y_bin, x_bin] += 1
        end

        # update progress meter
        if show_progress
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
        end
    end

    return phase_map_count

end



"""
    bin_2D_quantity_log!( phase_map_count, phase_map, 
                          x_q, y_q, bin_q, x_lim, y_lim; 
                          Nbins::Int=100,
                          show_progress::Bool=true)

Get a 2D histogram of `x_q` and `y_q` in the limits `x_lim`, `y_lim` over a number if bins `Nbins` in ``log``-space for pre-allocated arrays `phase_map_count` and `phase_map`.
Should be used if computing a 2D phase map over multiple files.
"""
function bin_2D_quantity_log!(phase_map_count, phase_map,
    x_q, y_q, bin_q, x_lim, y_lim;
    Nbins::Int=100, show_progress::Bool=true)

    # get logarithmic bin spacing
    dlogx = (log10(x_lim[2]) - log10(x_lim[1])) / Nbins
    dlogy = (log10(y_lim[2]) - log10(y_lim[1])) / Nbins

    # optional progress meter
    if show_progress
        P = Progress(size(x_q, 1))
        idx_p = 0
    end

    @inbounds for i = 1:size(x_q, 1)

        x_bin = 1 + floor(Int64, (log10(x_q[i]) - log10(x_lim[1])) / dlogx)
        y_bin = 1 + floor(Int64, (log10(y_q[i]) - log10(y_lim[1])) / dlogy)

        if (1 <= x_bin <= Nbins) && (1 <= y_bin <= Nbins)

            phase_map_count[y_bin, x_bin] += 1

            phase_map[y_bin, x_bin] += bin_q[i]

        end

        # update progress meter
        if show_progress
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
        end
    end

    return phase_map_count, phase_map

end

