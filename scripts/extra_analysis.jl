using AfricanMalariaVectors, Rasters, RasterDataSources
using DataFrames, Random, StatsBase
import SpeciesDistributionModels as SDM
using CairoMakie, Colors

const AMV = AfricanMalariaVectors

# Set a path for Raster data
# ENV["RASTERDATASOURCES_PATH"] = "my/data/path"

# These species names are used throughout
const FOCUS_SPECIES = (:gambiae, :arabiensis, :funestus, :moucheti_sl, :coluzzii, :nili_sl)

# This is where to write the actual results
global resultspath = "manuscript"

# These are all the bioclim variables we considered
considered_vars = ((Symbol("bio$i") for i in (1:7..., 10:17...))..., :gsl, :gsp, :gst, :npp)

# Prepare data - just current bioclim
bioclim = AMV.load_current_bioclim(considered_vars, false, 5)
land = AMV.get_land() # get polygons of countries within region of interest minus lakes
roi = AMV.get_roi(land, bioclim.bio12) # get region of interest (sub-saharan africa) as a raster
bioclim_msk = mask(bioclim; with = roi)

# load occurrence data
all_occurrences = load_occurrences()
occs_filtered = filter(x -> x.year_start < 2010 && x.year_end > 1980, all_occurrences)
occs_by_species = NamedTuple(
    S => filter(
        x -> Symbol(x.species) == S || (!ismissing(x.complex) && Symbol(x.complex) == S), occs_filtered
    ).geometry 
    for S in FOCUS_SPECIES
)
ocs_no_dupes = unique(occs_filtered, [:species, :geometry, :complex])

# Generate background points for each species
all_occs = extract(bioclim_msk, ocs_no_dupes.geometry, skipmissing = true, geometry=false)
occs_sample = sample(Xoshiro(1234), all_occs, 2000; replace = false)

# get random bg points
bio_sample = Rasters.sample(Xoshiro(1234),bioclim_msk, 1000; skipmissing = true, geometry = false)

bio_mat = hcat(
    reduce(hcat, collect.(bio_sample)),
    reduce(hcat, collect.(occs_sample))
)
bio_c = cor(bio_mat')

var_indices = findall(x -> x in (:bio1, :bio7, :bio12, :bio14, :bio15), considered_vars)

# Select the submatrix where both row and column are in selected_vars
bio_c_selected = bio_c[var_indices, var_indices]
findmax(abs, filter(x -> x != 1, bio_c_selected))

isselected = falses(19, 19)
for (i, v) in enumerate(considered_vars)
    if v in (:bio1, :bio7, :bio12, :bio14, :bio15)
        isselected[i,:] .= true
        isselected[:,i] .= true
    end
end

# Assume bio_c is your correlation matrix, and you have 19 variables
nvars = size(bio_c, 1)
varnames = collect(string.(considered_vars)) # or use actual variable names if available

corfig = let
    fig = Figure(size = (800, 700))
    ax = Axis(fig[1, 1], 
        xticks = (1:nvars, varnames), 
        yticks = (1:nvars, varnames), 
        xlabel = "Variable", ylabel = "Variable", 
        xticklabelrotation = pi/3
    )
    hm = heatmap!(ax, bio_c, colorrange = (-1, 1), colormap = :coolwarm)
    Colorbar(fig[1, 2], hm, label = "Correlation coefficient")
    heatmap!(ax, isselected; colormap = [RGBA(0,0,0,0), RGBA(0.5,0.5,0.5,0.5)], colorrange = (0,1))
    # Draw squares around selected cells
    for i in var_indices, j in var_indices
        i == j && continue # Skip diagonal
        y = [i-0.5, i-0.5, i+0.5, i+0.5, i-0.5]
        x = [j-0.5, j+0.5, j+0.5, j-0.5, j-0.5]
        lines!(ax, x, y, color=:black, linewidth=2)
    end

    # Overlay correlation values as text
    for i in 1:nvars, j in 1:nvars
        val = round(bio_c[i, j], digits=2)
        color = :black
    # color = abs(val) > 0.8 ? :red : :black  # Highlight high correlations
        text!(ax, j, i, text = string(val), color = color, fontsize = 14, align = (:center, :center))
    end
    fig
end

save(joinpath(resultspath, "bioclim_correlations.png"), corfig)

## ROI plot
nondesert = mask(bioclim.bio12; with = land) .> 200
roi_desert = (nondesert .* roi) .+ roi

roi_fig = let
    fig = Figure(size = (600, 600))
    ax = map_axis(fig[1,1])
    pl = plot!(
        fig[1, 1], roi_desert; 
        colormap = [:lightgray, :dimgray, :black], 
    )
    for x in occs_thinned
        scatter!(ax, x; color = :lightgray, markersize = 5)
    end
    legendelements = [
        PolyElement(color = :black),
        PolyElement(color = :dimgray),
        PolyElement(color = :lightgray),
        MarkerElement(marker = :circle, color = :lightgray, markersize = 5)
    ]
    legendlabels = ["Non-desert", "Buffer area", "Outside ROI", "Occurrence"]
    Legend(fig[1,1], legendelements, legendlabels; tellheight = false, tellwidth = false, halign = 0.15, valign = 0.2)
    fig
end

save(joinpath(resultspath, "roi_plot.png"), roi_fig)

### Number sorted out per species
using SummaryTables
import WriteDocx as W

by_species = NamedTuple(
    S => filter(
        x -> Symbol(x.species) == S || (!ismissing(x.complex) && Symbol(x.complex) == S), all_occurrences
    )
    for S in FOCUS_SPECIES
)

occ_counts = map(by_species,occs_thinned) do x, o
    post_1980 = filter(x -> x.year_end > 1980, x)
    pre_2010 = filter(x -> x.year_start < 2010, post_1980)   

    return (
        initial = nrow(x),    
        post_1980 = nrow(post_1980),
        pre_2010 = nrow(pre_2010),
        final = length(o)
    )
end

occ_counts_mat = reduce(hcat, collect.(values(occ_counts)))'
occ_counts_cells = Cell.(occ_counts_mat)
species_cells = Cell.("An. " .* as_label.(FOCUS_SPECIES), italic = true, halign = :left) |> collect
headers = Cell.(["Species", "Unfiltered", "Post-1980", "Pre-2010", "Spatial thinning"], bold = true, border_bottom = true)
table = SummaryTables.Table([headers'; species_cells occ_counts_cells])
doc = W.Document(
    W.Body([
        W.Section([
            W.Paragraph([
                W.Run([W.Text("Table S1")]),
            ]),
            SummaryTables.to_docx(table),
        ]),
    ]),
)
W.save(joinpath(resultspath, "occurrence_count_table.docx"), doc)

############################
#### Variable selection ####
############################
import Maxnet, MLJFlux, MLJModels, MLJBase
using SummaryTables, Printf
import WriteDocx as W

# Load lots of data
bioclim = AMV.load_current_bioclim(considered_vars, true, 5)
land = AMV.get_land() # get polygons of countries within region of interest minus lakes
roi = AMV.get_roi(land, bioclim.bio12) # get region of interest (sub-saharan africa) as a raster
lulc = AMV.load_lulc(roi).current

# combine bioclim and lulc into a single rasterstack
predictors = RasterStack(layers(bioclim)..., lulc)
# mask out the region of interest
predictors = mask(predictors; with = roi)

# extract occurrence data for each species
occs_bio = map(x -> extract(predictors, x; skipmissing = true, geometry = true), occs_by_species)
# thin occurrence data with a thinning distance of 10 km
occs_thinned = map(x -> SDM.thin(Xoshiro(0), x, 10_000), occs_bio)

# Generate background points for each species
all_occs = extract(predictors, ocs_no_dupes.geometry, skipmissing = true)
bgs = map(x -> generate_bg(Xoshiro(0), x, all_occs, predictors; fraction_uniform = 1, fraction_sampling = 1), occs_thinned)

# select models
models = (
    maxnet = Maxnet.MaxnetBinaryClassifier(features = "lqp"),
    neuralnetwork = MLJModels.Standardizer(count = true) |>
        MLJFlux.NeuralNetworkBinaryClassifier(builder = MLJFlux.MLP(hidden = (16, )), epochs = 50, batch_size = 20, rng = 0, embedding_dims = Dict(:lulc => 3)),
    gam = GAMClassifier(k = 5, gamma = 3)
)

resampled_data = SDM.sdmdata(x, bg; resampler, predictors = (:bio1, :bio12))
@profview ensemble = SDM.sdm(resampled_data, models, scitype_check_level = 0)

# evalaute models for each species
# use 5-fold grid resampling
resampler = AMV.GridCV(Xoshiro(0); nfolds = 5, sides = 3, geometry = predictors, strategy = :checkerboard)

cor_nt = NamedTuple{considered_vars}(NamedTuple{considered_vars}.(eachrow(bio_c)))

forward_selection_results = map(occs_thinned, bgs) do x, bg
    preds_set = (:lulc,)
    candidate_vars = considered_vars

    oldscore = 0.0
    while !isempty(candidate_vars)
        score, idx = findmax(candidate_vars) do newvar
            newset = (preds_set..., newvar, :geometry)
            resampled_data = SDM.sdmdata(x, bg; resampler, predictors = newset)
            ensemble = SDM.sdm(resampled_data, models, scitype_check_level = 0)
            ev = SDM.evaluate(ensemble; measures = (; AUC = SDM.auc))
            mean(ev.results.test.AUC.score)
        end
        score < oldscore && break
        oldscore = score
        pred_sel = candidate_vars[idx]
        preds_set = (preds_set..., pred_sel)
        candidate_vars = filter(x -> abs(cor_nt[pred_sel][x] < 0.7), candidate_vars)
    end
    return (predictors = preds_set, auc = oldscore)
end

sets = map(x -> x.predictors, forward_selection_results)
evs = map(sets, occs_thinned, bgs) do set, x, bg
    resampled_data = SDM.sdmdata(x, bg; resampler, predictors = (set..., :geometry))
    ensemble = SDM.sdm(resampled_data, models, scitype_check_level = 0)
    ev = SDM.evaluate(ensemble; measures = (; AUC = SDM.auc))
end

auc_cells = map(evs) do ev
    Cell(@sprintf "%0.2f (%0.2f)" mean(ev.results.test.AUC.score) mean(ev.results.train.AUC.score))
end |> collect
species_cells = Cell.("An. " .* as_label.(FOCUS_SPECIES), italic = true, halign = :left) |> collect
sets_cells = Cell.(join.(collect(sets), (", ",)); halign = :left)

headers = Cell.(["Species", "AUC", "Predictors"], bold = true, border_bottom = true)
table = SummaryTables.Table([headers'; species_cells auc_cells sets_cells])

doc = W.Document(
    W.Body([
        W.Section([
            W.Paragraph([
                W.Run([W.Text("Table S2")]),
            ]),
            SummaryTables.to_docx(table),
        ]),
    ]),
)
W.save(joinpath(resultspath, "forward_selection_results.docx"), doc)
