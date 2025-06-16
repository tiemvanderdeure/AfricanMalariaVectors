using AfricanMalariaVectors
using Rasters, RasterDataSources
import Maxnet, MLJFlux, MLJModels, MLJBase
using Random, StatsBase, Tables, CSV, Dates, DataFrames
import SpeciesDistributionModels as SDM
using ProgressMeter, ThreadsX

const AMV = AfricanMalariaVectors

# Set a path for Raster data
# ENV["RASTERDATASOURCES_PATH"] = "my/data/path"

# These species names are used throughout
const FOCUS_SPECIES = (:gambiae, :arabiensis, :funestus, :moucheti_sl, :coluzzii, :nili_sl)

# This is where to write the actual results
global resultspath = "manuscript"

# Load lots of data
bioclim = load_bioclim((:bio1, :bio7, :bio12, :bio14, :bio15))
land = AMV.get_land() # get polygons of countries within region of interest minus lakes
roi = AMV.get_roi(land, bioclim.current.bio12) # get region of interest (sub-saharan africa) as a raster
lulc = AMV.load_lulc(roi)
pop = load_worldpop(roi)

# combine bioclim and lulc into a single rasterstack
predictors = map((b, l) -> RasterStack(layers(b)..., l), bioclim, lulc)
# mask out the region of interest
predictors = map(b -> mask(b; with = roi), predictors)

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

# extract occurrence data for each species
occs_bio = map(x -> extract(predictors.current, x; skipmissing = true, geometry = true), occs_by_species)
# thin occurrence data with a thinning distance of 10 km
occs_thinned = map(x -> SDM.thin(Xoshiro(0), x, 10_000), occs_bio)

# Generate background points for each species
all_occs = extract(predictors.current, ocs_no_dupes.geometry, skipmissing = true)
bgs = map(x -> generate_bg(Xoshiro(0), x, all_occs, predictors.current; fraction_uniform = 1, fraction_sampling = 1), occs_thinned)

# select models
models = (
    maxnet = Maxnet.MaxnetBinaryClassifier(features = "lqp"),
    neuralnetwork = MLJModels.Standardizer(count = true) |>
        MLJFlux.NeuralNetworkBinaryClassifier(builder = MLJFlux.MLP(hidden = (16, )), epochs = 50, batch_size = 20, rng = 0, embedding_dims = Dict(:lulc => 3)),
    gam = GAMClassifier(k = 5, gamma = 3)
)

# evalaute models for each species
# use 5-fold grid resampling
resampler = AMV.GridCV(Xoshiro(0); nfolds = 5, sides = 3, geometry = predictors.current, strategy = :checkerboard)
measures = (
    AUC = SDM.auc, 
    TSS = SDM.StatisticalMeasures.BalancedAccuracy(true),
    CBI = AMV.cbi(), 
)

# fit and evaluate for each species
evaluations = @showprogress map(occs_thinned, bgs) do x, bg
    resampled_data = SDM.sdmdata(x, bg; resampler)
    ensemble = SDM.sdm(resampled_data, models, scitype_check_level = 0)
    SDM.evaluate(ensemble; measures)
end;

## fit the models and construct ensembles
ensembles = @showprogress map(occs_thinned, bgs) do x, bg
    data = SDM.sdmdata(x, bg)
    SDM.sdm(data, models, scitype_check_level = 0)
end;

## generate predictions for current and future climate
preds = map(predictors) do b
    ThreadsX.map(ensembles) do e
        SDM.predict(e, b; reducer = mean)
    end |> NamedTuple{keys(ensembles)}
end

preds_future_mean = map(preds.future) do p
    dropdims(mean(p; dims = :gcm); dims = :gcm)
end

## thresholds for each species, using 10th percentile presence
thresholds = map(ensembles, occs_thinned) do e, x
    quantile(SDM.predict(e, x; reducer = mean), 0.1)
end

# make binary maps for each species
preds_binary = map(p -> map((p, t) -> p .> t, p, thresholds), preds)

## population at risk
pop_at_risk = (
    current = map(p -> sum(skipmissing(pop .* p)), preds_binary.current), 
    future = map(preds_binary.future) do p
            [sum(skipmissing(pop .* s)) for s in eachslice(p; dims = (:gcm, :ssp, :Ti))]
        end
)

#### Malaria data
malariadata = AMV.load_malaria()

# extract suitabilities for each vector species at each survey location
for S in keys(preds.current)
    malariadata[!, S] =
        getindex.(
            extract(preds.current[S], malariadata, geometrycolumn = (:Long, :Lat), skipmissing = false, geometry = false),
            1
        )
end
dropmissing!(malariadata, collect(keys(preds.current)))

grp = groupby(malariadata, collect(keys(preds.current)))
malariadata_by_pixel = DataFrames.combine(grp, :PfPR2_10 => mean => :PfPR2_10, :Lat => first => :Lat, :Long => first => :Long)
malaria_correlations = map(FOCUS_SPECIES) do S
    cor(malariadata_by_pixel[!, :PfPR2_10], malariadata_by_pixel[!, S])
end |> NamedTuple{FOCUS_SPECIES}

#### Model explanations
# use shapley values to explain the ensembles
using Shapley
# Takes 15 minutus or so
shapvals = map(ensembles) do e
    AMV.shapley(
        Shapley.MonteCarlo(Shapley.CPUThreads(), 32), 
        e,
        predictors.current
    )
end

variable_importances = map(s -> maplayers(x -> mean(abs, skipmissing(x)), s), shapvals)

# this takes ~ 15 minutes
bshapvals_future = ThreadsX.map(ensembles) do e
    rast = cat(
        (
            Bshap(x -> SDM.predict(e, x; reducer = mean), predictors.current, fut) for 
            fut in eachslice(predictors.future[ssp = At(SSP370), Ti = At(Date(2085))]; dims = :gcm)
        )...;
        dims = dims(predictors.future, :gcm)
    )
    dropdims(mean(rast; dims = :gcm); dims =:gcm)
end |> NamedTuple{keys(ensembles)}


# These numbers are cited in the text
open(joinpath(resultspath, "in_text_numbers.txt"), "w") do io
    gamb_end_of_cent = pop_at_risk.future.gambiae[Ti = 2, ssp = 2] ./ 1e6
    nili_end_of_cent = pop_at_risk.future.nili_sl[Ti = 2, ssp = 2] ./ 1e6
    println(io, "Nili sl end of century:  $(mean(nili_end_of_cent)), $(extrema(nili_end_of_cent))")
    println(io, "Gambiae end of century:  $(mean(gamb_end_of_cent)), $(extrema(gamb_end_of_cent))")
end