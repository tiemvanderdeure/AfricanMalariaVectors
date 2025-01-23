using CairoMakie#, GeoMakie
using Dates, Printf, Tables
using SpeciesDistributionModels: Loess
using NaturalEarth
countries = naturalearth("ne_10m_admin_0_countries")

global imagepath = joinpath("images", "nat_cc")

#set_theme!(Theme(font = :))
#ytheme = Theme(font = "Helvetica")

## New Figure 1

## Figure with scatter plots of temperature and precipitation
begin
    temp_prec = bioclim.current[(:bio1, :bio12)]
    temp_prec_points = map(x -> (extract(temp_prec, x; geometry = false, skipmissing = true)), occs_thinned)
    temp_prec_points_all = vcat(temp_prec_points...) |> Tables.columntable
    temp_prec_focus_species = map(Tables.columntable, temp_prec_points[FOCUS_SPECIES])
    points_all = vcat(map(p -> getfield.(p, :geometry), occs_thinned)...)
    points_focus_species = map(p -> getfield.(p, :geometry), occs_thinned[FOCUS_SPECIES])

    fig = Figure(size = (1640, 1500), fontsize = 16)
    for (i, S) in enumerate(FOCUS_SPECIES)
        gl = fig[mod1(i, 3), (i-1) ÷ 3] = GridLayout()

        data = (temp_prec_points_all, temp_prec_focus_species[S])
        pointdata = (points_all, points_focus_species[S])
        colors = (:lightblue, :purple)
        labels = ("All\nAnopheles", "A. $(as_label(S))")

        ## map on the left
        ma = map_axis(gl[1,1], tellwidth = true, tellheight = true, aspect = nothing)
        poly!(ma, countries.geometry; color = :transparent, strokewidth = 0.3, strokecolor = :grey)
        for (p, col) in zip(pointdata, colors)
            scatter!(ma, p, markersize = 3, color = col)
        end
        # force map element to be square
        colsize!(gl, 1, Aspect(1, 1))

        ## scatter and density on the right
        glright = gl[1, 2] = GridLayout()
        axmain = Axis(glright[2,1]; xlabel = "Annual Mean Temperature (°C)", ylabel = "Annual Precipitation (mm)")
        axtop = Axis(glright[1,1], limits = (13, 32, 0, nothing), ygridvisible = false)
        axright = Axis(glright[2,2], limits = (0, nothing, 0, 4000), xgridvisible = false)
        linkyaxes!(axmain, axright)
        linkxaxes!(axmain, axtop)

        for (d, col, lab) in zip(data, colors, labels)
            scatter!(axmain, d.bio1, d.bio12; color = col, markersize = 3, label = lab)
            density!(axtop, d.bio1; color = (col, 0.7))
            density!(axright, d.bio12 ; color = (col, 0.7), direction = :y)
        end
        #leg = Legend(glright[1, 2], axmain, tellheight = true)

        # layout tweaks
        rowsize!(glright, 1, Relative(0.15)) # upper density plot fills 15%
        colsize!(glright, 2, Aspect(1, 1)) # right density plot is as wide as upper is high
        colsize!(glright, 1, Aspect(2, 1)) # scatter plot is square
        colgap!(glright, 5)
        rowgap!(glright, 5)
        hidedecorations!(axtop, grid = false)
        hidedecorations!(axright, grid = false)

        # label on top
        Label(gl[-1, 1:2], "Anopheles " * as_label(S), font = :bold_italic, fontsize = 20)
        Label(gl[0, 1:2], "n = $(length(pointdata[2]))")
        letter_label(gl[-1:0,1:2], i; halign = :left)
        rowgap!(gl, 2)
    end

    [rowsize!(fig.layout, i, Fixed(400)) for i in 1:3] # fix the rowsize so layout can be resolved
    resize_to_layout!(fig)
    save(joinpath(imagepath, "figure1.png"), fig)
end

#### Figure 2 - maps with projections
let cmap = Reverse(:Spectral), colorrange = (0,1)
    for ssp in SSPS
        fig = Figure()
        for (i, s) in enumerate(FOCUS_SPECIES)
            axes = [map_axis(fig[j, i]; land) for j in 1:3]
            plot!(axes[1], preds.current[s]; colormap = cmap, colorrange)
            plot!(axes[2], preds_future_mean[s][date = At(Date(2055)), ssp = At(ssp)]; colormap = cmap, colorrange)
            plot!(axes[3], preds_future_mean[s][date = At(Date(2085)), ssp = At(ssp)]; colormap = cmap, colorrange)
            Label(fig[1, i, Top()], "A. " * as_label(s), font = :bold_italic, padding = (0,0,5,0))
            colsize!(fig.layout, i, Aspect(1, 1.0))
        end
        for (i, l) in enumerate(["current", "2041-2070", "2071-2100"])
            Label(fig[i,0], l; tellwidth = true, tellheight = false, font = :bold, rotation = π/2)
        end
        Colorbar(fig[1:3,length(FOCUS_SPECIES)+1]; label = "Climatic suitability", colormap = cmap, colorrange)
        rowgap!(fig.layout, 5); colgap!(fig.layout, 5)
        resize_to_layout!(fig)
        filename = ssp == SSP370 ? "figure2.png" : "extended_figure1.png"
        save(joinpath(imagepath, filename), fig)
    end
end

#### Figure 3 - Shapley values
# data - combine temperature and precipitation variables
shapvals_comb = map(shapvals) do s
    RasterStack(
        (temperature = s.bio1 .+ s.bio7,
        precipitation = s.bio12 .+ s.bio14 .+ s.bio15,
        lulc = s.lulc))
end    
bshapvals_future_comb = map(bshapvals_future) do s
    RasterStack(
        (temperature = s.bio1 .+ s.bio7,
        precipitation = s.bio12 .+ s.bio14 .+ s.bio15,
        lulc = s.lulc))
end

bg = Rasters.sample(Xoshiro(0), bioclim.current, 10000; skipmissing = true, geometry = (X, Y)) |> Tables.columntable
bg_combined = (temperature = bg.bio1, precipitation = bg.bio12 , lulc = bg.lulc)

# labels and colors
xlabels = (
    temperature = "Annual Mean Temperature (°C)",
    precipitation = "Annual Precipitation (mm)",
    lulc = "",
)
xlimits = (
    temperature = (nothing, nothing),
    precipitation = (0, 3800),
    lulc = (nothing, nothing)
)
titles = (
    temperature = "Temperature",
    precipitation = "Precipitation",
    lulc = "Land Use"
)

let colorrange = (-0.5, 0.5),
    colormap = :vik,
    xlabels = xlabels, xlimits = xlimits, titles = titles,
    shapvals_comb = shapvals_comb, bshapvals_future_comb = bshapvals_future_comb,
    bg = bg, bg_combined = bg_combined

    for (j, K) in enumerate(keys(shapvals))
        fig = Figure(size = (1000, 750), alignmode = Mixed(left = 0))

        # get shapley values at 10_000 random points for the scatter plot
        bg_shapvals = extract(shapvals_comb[K], bg, skipmissing = true) |> Tables.columntable

        for (i, l) in enumerate(layers(shapvals_comb[K]))
            var = name(l)
            Label(fig[0, i], titles[var], font = :bold, tellwidth = false)

            # map
            ax = map_axis(fig[1, i])
            plot!(ax, l; colorrange, colormap)

            # scatter
            ax2 = Axis(fig[2, i], xlabel = xlabels[var], limits = (xlimits[var]..., nothing, nothing))
            if i == 1
                ax2.ylabel = "Contribution to suitability \n (Shapley value)"
            else
                linkyaxes!(ax2, content(fig[2,1]))
                ax2.yticklabelsvisible = false
                ax2.yticksvisible = false
            end
            if var == :lulc
                ax2.xticks = (1:length(levels(bg_combined[var])), levels(bg_combined[var]))
                ax2.xticklabelrotation = pi/6
                boxplot!(ax2, bg_combined[var].refs, bg_shapvals[var]; color = :grey, markersize = 2, whiskerlinewidth = 3)
            else
                scatter!(ax2, bg_combined[var], bg_shapvals[var]; color = :black, markersize = 1)
            end

            # map future
            ax = map_axis(fig[3, i])
            plot!(ax, bshapvals_future_comb[K][var]; colorrange, colormap)
        end
        # add colorbars
        Colorbar(fig[1,1, Left()]; label = "Contribution to suitability", colormap, colorrange, flipaxis =false, halign = 0.9,
            height = Relative(0.7))
        Colorbar(fig[3,1, Left()]; label = "Contribution to suitability change", colormap, colorrange, flipaxis =false, halign = 0.9,
            height = Relative(0.7))
        # force maps to be aspect 1
        for i in 1:3 colsize!(fig.layout, i, Aspect(1,1)) end
        colgap!(fig.layout, 10)
        resize_to_layout!(fig)

        # save figures
        if K == :gambiae
            save(joinpath(imagepath, "figure3.png"), fig)
        else
            save(joinpath(imagepath, "extended_figure_$(j)_$K.png"), fig)
        end
    end
end


#### Figure 4 - Malaria records and correlations
fig4 = begin
    fig = Figure(size = (1200, 800), fontsize = 16)
    glleft = fig[1, 1] = GridLayout(alignmode = Outside())
    glright = fig[1, 2] = GridLayout(alignmode = Outside())

    mapaxis = map_axis(glleft[1,1])
    poly!(mapaxis, countries.geometry, color = :transparent, strokecolor = :black, strokewidth = 0.3)
    sc = scatter!(
        malariadataclean.Long, malariadataclean.Lat; 
        color = malariadataclean.PfPR2_10, colorrange = (0,100), colormap = :viridis,
        markersize = 5, strokewidth = 0
    )
    Colorbar(
        glleft[1,1], sc, 
        valign = 0.15, halign = 0.2, height = Relative(0.4), tellwidth = false, tellheight = false,
        label = rich(rich("P. falciparum", font = :italic), " rate")
    )

    for (i, S) in enumerate(FOCUS_SPECIES)
        r = mod1(i, 3) # row index
        c = (i-1) ÷ 3 + 1 # col index
        ax = Axis(
            glright[r, c]; 
            title = "Anopheles " * as_label(S), titlefont = :bold_italic,
            xlabel = "Suitability", ylabel = rich(rich("P. falciparum", font = :italic), " rate"),
            xticklabelsvisible = r == 3, yticklabelsvisible = c == 1,
            xlabelvisible = r == 3, ylabelvisible = c == 1,
            xticksvisible = r == 3, yticksvisible = c == 1,
            )

        if i > 1
            linkyaxes!(ax, content(glright[1,1]))
            linkyaxes!(ax, content(glright[1,1]))
        end

        windowwidth = 0.2
        start = 0.0
        step = 0.05
        stop = 1.0 - windowwidth
        binned_values = map(start:step:stop) do i
            malariadataclean.PfPR2_10[i .<= malariadataclean[!, S] .< (i + windowwidth)]
        end

        xs = ((start:step:stop) .+ windowwidth / 2)[.~isempty.(binned_values)]
        binned_values = filter(!isempty, binned_values)

        quantiles = quantile.(binned_values, Ref([0.25, 0.75]))
        medians = median.(binned_values)

        size = sqrt.(length.(binned_values))
        scatter!(ax, xs, medians; color = :black, markersize = size / 6)
        errorbars!(ax, xs, medians, medians .- first.(quantiles), last.(quantiles) .- medians; color = :black, linewidth = size / 20)

        text!(ax, 0.1, 60, text = (@sprintf "cor: %.2f" malaria_correlations[S]))
    end

    # Add a legend for markersize
    ns = [10, 50, 100, 250]
    group_size = [MarkerElement(marker = :circle, color = :black,
        strokecolor = :transparent,
        markersize = sqrt(n)) for n in ns]
    Legend(
        glleft[1,1], group_size, string.(ns), "Surveys", 
        tellheight = false, tellwidth = false,
        halign = 1.13, valign = :bottom,
        framevisible = true, rowgap = 0, titlegap = 0
    )

    # add a and b labels
    for i in 1:2 letter_label(fig[1,i], i; halign = :left, valign = :top, fontsize = 22) end

    # layout tweaks
    rowgap!(glright, 10); colgap!(glright, 10)
    rowsize!(fig.layout, 1, Aspect(1,1))
    colgap!(fig.layout, 10)
    resize_to_layout!(fig)

    save(joinpath(imagepath, "figure4.png"), fig)
    fig
end

#### Figure 5 - population and suitability maps
# get a bivariate color map
xticks = [0, 0.25, 0.5, 0.75, 1]
yticks = [0, 50, 500, 5000, Inf]
brewer_seqseq2 = [
    colorant"#ffac36" colorant"black"
    colorant"#f3f3f3" colorant"#209ebe"
]

#=
tolochko_redblue = [
    colorant"#dd0027" colorant"#4f2d4c"
    colorant"#dcdcdc" colorant"#0072a9"
]
tolochko_redblue2 = [
    colorant"#4f2d4c" colorant"#dd0027"
    colorant"#dcdcdc" colorant"#0072a9"
]
=#

fig5 = let bivcmap = bivariate_colormap(xticks, yticks; colors = brewer_seqseq2), 
    barscmap = :blues,
    s = :gambiae

    fig = Figure(size = (600, 800), fontsize = 12)

    mapgl = fig[1,1:2] = GridLayout(alignmode = Mixed(left = 0, top = 0), valign = :top)
    for (j, ssp) in enumerate(SSPS)
        ax = map_axis(mapgl[j, 1], title = j == 1 ? "Current" : "")
        plot!(ax, rasters_to_cmap(preds.current[s], pop; cmap = bivcmap))
        for (i, date) in enumerate(DATES)
            ax = map_axis(mapgl[j, i+1], title = j == 1 ? as_label(date) : "")
            plot!(ax, rasters_to_cmap(preds_future_mean[s][date = At(date), ssp = At(ssp)], pop; cmap = bivcmap))
        end
        Label(mapgl[j, 1, Left()], string(ssp), rotation = pi/2, tellheight = false, font = :bold)
        rowsize!(mapgl, j, Aspect(1, 1))
    end
    # tweaks
    colgap!(mapgl, 5)
    rowgap!(mapgl, 5)

    ## Bars
    bargl = fig[2,1] = GridLayout(alignmode = Mixed(left = 0, top = 0))

    for (i, ssp) in enumerate(SSPS)
        ax = Axis(bargl[1,i], 
            limits = (nothing, nothing, 0, nothing),
            xticks = (eachindex(FOCUS_SPECIES), collect(as_label.(FOCUS_SPECIES))), xticklabelrotation = pi/4, 
            yticksvisible = i == 1, yticklabelsvisible = i == 1,
            ylabel = i == 1 ? "population at risk (millions)" : "",
            ytickformat = x -> string.(floor.(Int, x ./ 1_000_000)),
            xgridvisible = false,
            title = string(ssp)
        )
        hidespines!(ax, :t, :r)
        if i > 1
            linkaxes!(ax, content(bargl[1,1]))
        end
        data = vcat(([pop_at_risk.current[s]; vec(mean(pop_at_risk.future[s][ssp = At(ssp)]; dims = :gcm))] for s in FOCUS_SPECIES)...)
        cats = repeat(eachindex(FOCUS_SPECIES), inner = 3)
        dodge = repeat(1:3, outer = length(FOCUS_SPECIES))
        barplot!(ax, cats, data; dodge, color = dodge, dodge_gap = 0, colormap = barscmap)
        # add uncertainty bars
        data_whiskers = vcat(([(); vec(extrema(pop_at_risk.future[s][ssp = At(ssp)]; dims = :gcm))] for s in FOCUS_SPECIES)...)
        AMV.myerrorbars!(ax, cats, data_whiskers; dodge, dodge_gap = 0, color = :black, whiskerwidth = 5, linewidth = 1)
    end
    colgap!(bargl, 5)

    ## Legends in the bottom right
    legendgl = fig[2,2] = GridLayout()
    # maps legend
    bivlegendlayout = legendgl[1,1] = GridLayout(
        halign = :right, valign = 1.1, tellwidth = false, tellheight = false, alignmode = Outside(), width = Relative(0.95)
    )
    rowsize!(bivlegendlayout, 1, Aspect(1, 1.0))

    legendpopf(ticks)= [string.(Int.(ticks[1:end-1])); ">"* string(Int(ticks[end-1]))]
    bivariate_cmap_legend(
        bivlegendlayout[1,1], bivcmap; 
        yticksf = legendpopf, ylabel = "Population\ndensity", xlabel = "A. gambiae \nsuitability",
        xlabelpadding = 0, ylabelpadding = 0)

    # bars legend
    labels = ["current", as_label.(DATES)...]
    elements = [PolyElement(color = i, colorrange = (1, length(labels)), colormap = barscmap) for i in eachindex(labels)]

    Legend(fig[2,2], elements, labels, halign = :left, valign = :bottom, framevisible = false,
        padding = (0,0,0,0))

    # labels
    for i in 1:2 letter_label(fig[i,1], i; halign = :left, valign = :top, fontsize = 16) end

    ## global tweaks
    colsize!(fig.layout, 2, Relative(0.28))
    rowsize!(fig.layout, 2, Fixed(200))
    colsize!(fig.layout, 1, Aspect(2, 2.0))
    rowsize!(fig.layout, 1, Aspect(1, 1.0))
    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 5)
    resize_to_layout!(fig)
    save(joinpath(imagepath, "figure5.png"), fig)
    fig
end

# These numbers are cited in the text
gamb_end_of_cent = pop_at_risk.future.gambiae[date = 2, ssp = 2]
nili_end_of_cent = pop_at_risk.future.nili_sl[date = 2, ssp = 2]
mean(nili_end_of_cent)
extrema(nili_end_of_cent)
mean(gamb_end_of_cent)
extrema(gamb_end_of_cent)