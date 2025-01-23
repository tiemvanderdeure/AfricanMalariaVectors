using AfricanMalariaVectors
using Rasters, RasterDataSources
import Maxnet, MLJFlux, MLJModels, MLJBase
using Random, StatsBase, Tables, CSV, Dates, DataFrames
import SpeciesDistributionModels as SDM
using ProgressMeter, ThreadsX

const AMV = AfricanMalariaVectors

FOCUS_SPECIES = (:gambiae, :arabiensis, :funestus, :moucheti_sl, :coluzzii, :nili_sl)

PREDICTORS = (:bio1, :bio7, :bio12, :bio14, :bio15)
bioclim = load_bioclim(PREDICTORS)

land = AMV.get_land() # get polygons of countries within region of interest minus lakes
roi = AMV.get_roi(land, bioclim.current.bio12) # get region of interest (sub-saharan africa) as a raster
lulc = AMV.load_lulc(roi)
pop = load_worldpop(roi)

bioclim = map((b, l) -> RasterStack(layers(b)..., l), bioclim, lulc)
bioclim = map(b -> mask(b; with = roi), bioclim)

# load occurrence data
all_occurrences = load_occurrences()

occs_by_species = NamedTuple(
    S => filter(
        x -> Symbol(x.species) == S || (!ismissing(x.complex) && Symbol(x.complex) == S), all_occurrences
    ).geometry 
    for S in FOCUS_SPECIES
)

ocs_no_dupes = unique(all_occurrences, [:species, :geometry, :complex])

# extract occurrence data for each species
occs_bio = map(x -> extract(bioclim.current, x; skipmissing = true, geometry = true), occs_by_species)
# thin occurrence data with a thinning distance of 10 km
occs_thinned = map(x -> SDM.thin(x, 10_000), occs_bio)

# Generate background points for each species
all_occs = extract(bioclim.current, ocs_no_dupes.geometry, skipmissing = true)
bgs = map(x -> generate_bg(Xoshiro(0), x, all_occs, bioclim.current; fraction_uniform = 1, fraction_sampling = 1), occs_thinned)

# select models
models = (
    maxnet = Maxnet.MaxnetBinaryClassifier(features = "lqp"),
    neuralnetwork = MLJModels.Standardizer(count = true) |> #MLJModels.OneHotEncoder() |> # much faster than without!
        MLJFlux.NeuralNetworkBinaryClassifier(builder = MLJFlux.MLP(hidden = (16, )), epochs = 50, batch_size = 20, rng = 0, embedding_dims = Dict(:lulc => 3)),
    gam = GAMClassifier(k = 5, gamma = 3)
)

# evalaute models for each species
# use 4-fold grid resampling
resampler = AMV.GridCV(Xoshiro(0); nfolds = 5, sides = 3, geometry = bioclim.current, strategy = :checkerboard)
measures = (
    AUC = SDM.auc, 
    TSS = SDM.StatisticalMeasures.BalancedAccuracy(true),
    CBI = AMV.cbi(), 
)

#x = occs_thinned.moucheti_sl; bg = bgs.moucheti_sl
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
preds = map(bioclim) do b
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
            [sum(skipmissing(pop .* s)) for s in eachslice(p; dims = (:gcm, :ssp, :date))]
        end
)

#### Malaria data
malariadata = AMV.load_malaria()
rename!(malariadata, Symbol("PfPR2-10") => :PfPR2_10)
malariapointdata = filter(
    r -> r.AREA_TYPE == "Point" && r.YY >= 1980 && r.YY <= 2010, 
    malariadata
)

malariadataclean = select(malariapointdata, [:Lat, :Long, :PfPR2_10])
# extract suitabilities for each vector species at each survey location
for S in keys(preds.current)
    malariadataclean[!, S] =
        getindex.(
            extract(preds.current[S], malariapointdata, geometrycolumn = (:Long, :Lat), skipmissing = false, geometry = false),
            1
        )
end
dropmissing!(malariadataclean, collect(keys(preds.current)))

grp = groupby(malariadataclean, collect(keys(preds.current)))
malariadata_by_pixel = DataFrames.combine(grp, :PfPR2_10 => mean => :PfPR2_10, :Lat => first => :Lat, :Long => first => :Long)
malaria_correlations = map(FOCUS_SPECIES) do S
    cor(malariadata_by_pixel[!, :PfPR2_10], malariadata_by_pixel[!, S])
end |> NamedTuple{FOCUS_SPECIES}

#### Model explanations
# use shapley values to explain the ensembles
using Shapley
shapvals = map(ensembles) do e
    AnophelesSDMs.shapley(
        Shapley.MonteCarlo(Shapley.CPUThreads(), 32), 
        e,
        bioclim.current
    )
end

variable_importances = map(s -> map(x -> mean(abs, skipmissing(x)), s), shapvals)

# this takes ~ 15 minutes
bshapvals_future = ThreadsX.map(ensembles) do e
    rast = cat(
        (
            Bshap(x -> SDM.predict(e, x; reducer = mean), bioclim.current, fut) for 
            fut in eachslice(bioclim.future[ssp = At(SSP370), date = At(Date(2085))]; dims = :gcm)
        )...;
        dims = dims(bioclim.future, :gcm)
    )
    dropdims(mean(rast; dims = :gcm); dims =:gcm)

end |> NamedTuple{keys(ensembles)}
