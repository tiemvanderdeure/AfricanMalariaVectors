using Makie, Colors
using Rasters.Lookups

function map_axis(gp, args...; land = nothing, kw...)
    ax = Axis(gp, args...; limits = Tuple(LIMITS), aspect = 1, kw...)  
    hidedecorations!(ax); hidespines!(ax)
    return ax
end

# bivariate color map utils
function rasters_to_cmap(x, y; cmap)
    map(x, y) do xval, yval
        if ismissing(xval) || ismissing(yval) || xval === missingval(x) || yval === missingval(y) 
            RGBA(0,0,0,0) # transparent
        else
            cmap[X = Contains(xval), Y = Contains(yval)]
        end
    end
end

function bivariate_cmap_legend(
    gp, cmap; 
    halign = :center, valign = :center, 
    height = nothing, width = nothing, 
    xticksf = Base.Fix1(Broadcast.broadcast, string), yticksf = Base.Fix1(Broadcast.broadcast, string),
    kw... # to pass to axis
)
    xspan, yspan = val.(Rasters.span(dims(cmap)))
    xticks = [xspan[1]; xspan[2,:]]
    yticks = [yspan[1]; yspan[2,:]]
    ax = Axis(
        gp; 
        halign, valign, height, width, 
        tellheight = false, tellwidth = false,
        xticks = (0:size(xspan,2), xticksf(xticks)),
        yticks = (0:size(yspan,2), yticksf(yticks)),
        kw...
    )
    heatmap!(ax, 0:size(xspan,2), 0:size(yspan,2), parent(cmap))
end

function bivariate_colormap(xticks, yticks; colors)
    palette = colors2d.(range(0,1, length = length(xticks)-1), collect(range(0,1, length = length(yticks)-1))'; colors)
    ds = map((xticks, yticks), (X,Y)) do ticks, D
        D(Sampled(ticks[1:end-1]; span = Explicit(rotl90([ticks[2:end] ticks[1:end-1]]))))
    end
    DimArray(palette, ds)
end

function colors2d(x, y; colors)
    x1 = colors[2,2] * x + colors[2,1] * (1-x)
    x2 = colors[1,2] * x + colors[1,1] * (1-x)
    return x2 * y + x1 * (1-y)
end

# to add errorbars to a bar plot whiskers
function myerrorbars!(ax, x, y; dodge, dodge_gap = Makie.Automatic(), n_dodge = Makie.Automatic(), gap = 0.2, width = Makie.Automatic(), kw...)
    x̂, barwidth = Makie.compute_x_and_width(x, width, gap, dodge, n_dodge, dodge_gap)
    mask = .!isempty.(y)
    y = y[mask]
    lower = first.(y)
    upper = last.(y)
    errorbars!(ax, x̂[mask], first.(y), zeros(length(lower)), upper .- lower; kw...)
end

function compute_x(x, width, gap, dodge, dodge_gap)
    scale_width(dodge_gap, n_dodge) = (1 - (n_dodge - 1) * dodge_gap) / n_dodge
    function shift_dodge(i, dodge_width, dodge_gap)
        (dodge_width - 1) / 2 + (i - 1) * (dodge_width + dodge_gap)
    end
    width *= 1 - gap
    n_dodge = maximum(dodge)
    dodge_width = scale_width(dodge_gap, n_dodge)
    shifts = shift_dodge.(dodge, dodge_width, dodge_gap)
    return x .+ width .* shifts
end


# utils
as_label(S) = replace(string(S), "_sl" => " s.l.")
as_label(d::Date) = d == Date(2055) ?  "2041-2070" : "2071-2100"

letter_label(gp::GridPosition, i::Integer; fontsize = 20, font = :bold, tellheight = false, tellwidth = false, kw...) = 
    Label(gp, string('a'+i-1); font, fontsize, tellheight, tellwidth, kw...)



