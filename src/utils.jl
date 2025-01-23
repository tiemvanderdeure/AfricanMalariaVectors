using Random

sample_rstack(rs::AbstractRasterStack, i::Integer; kw...) = sample_rstack(Random.GLOBAL_RNG, rs, i; kw...)

function sample_rstack(
    rng::Random.AbstractRNG, rs::AbstractRasterStack, i::Integer; 
    wts = ones(dims(rs)), skipmissing::Bool = false, weigh_cellsize::Bool = false,
    replace = true, ordered = false, geometry = false)

    if skipmissing
        # generate a missing mask
        missingmask = boolmask(rs)
        wts .*= missingmask
    end

    if weigh_cellsize
        wts .*= Rasters.cellsize(rs)
    end

    weights = StatsBase.weights(wts)
    dimindices = Rasters.DimIndices(rs)
    indices = sample(rng, 1:prod(size(dimindices)), weights, i; replace = replace, ordered = ordered)
    v = map(idx -> getindex(rs, dimindices[idx]...), indices)

    if geometry
        dimpoints = DimPoints(rs)
        points = map(idx -> getindex(dimpoints, idx), indices)
        v = map(v, points) do v, point
            merge((; geometry = point), v)
        end
    end

    return v
end
#=
import Maxnet
maxn(args...; kw...) = Maxnet.MaxnetBinaryClassifier(args...; kw...)

import MLJFlux
ann(args...; kw...) = MLJFlux.NeuralNetworkClassifier(args...; kw...)
=#

## Spatial resampling
import Rasters.Extents
struct GridCV <: MLJBase.ResamplingStrategy
    nfolds::Integer
    samplinggrid::Raster{<:Integer}
end

GridCV(; kw...) = GridCV(Random.GLOBAL_RNG; kw...)
GridCV(rng::Random.AbstractRNG; geometry, nfolds = 4, sides = 5, strategy = :rand) = GridCV(rng, nfolds, sides, Extents.extent(geometry), strategy)
function GridCV(rng::Random.AbstractRNG, nfolds::Integer, sides::Integer, extent::Extents.Extent, strategy::Symbol)
    xbounds = range(extent.X...; step = sides)
    ybounds = range(extent.Y...; step = sides)
    resamplinggrid = Raster(
        ones(typeof(nfolds), 
            X(xbounds; sampling = DD.Intervals(DD.Start()), span = DD.Regular(sides)), 
            Y(ybounds; sampling = DD.Intervals(DD.Start()), span = DD.Regular(sides))
        );
        name = :fold
    )

    resamplinggrid .= eval(strategy)(rng, Base.OneTo(nfolds), size(resamplinggrid)...)
    return GridCV(nfolds, resamplinggrid)
end

function checkerboard(rng, values, nrows, ncols)
    nfolds = length(values)
    mod1.((1:nrows) .+ (1:ncols)', nfolds)
end

function MLJBase.train_test_pairs(resampler::GridCV, rows, X, y)
    fold = getfield.(Rasters.extract(resampler.samplinggrid, X), :fold)
    nfolds = maximum(fold) #resampler.nfolds
    map(1:nfolds) do f
        test = rows[fold .== f]
        train = Base.setdiff(rows, test)
        return (train, test)
    end
end

function generate_bg(rng, x, all_occs, bioc, bm = Rasters.boolmask(bioc); fraction_uniform, fraction_sampling)
    n_p = length(x)
    other_occs = Base.setdiff(all_occs, x)
    uni_bg = sample_rstack(rng, bioc, floor(Int, n_p * fraction_sampling); wts = bm, geometry = true)
    sampl_bg = sample(rng, other_occs, floor(Int, n_p * fraction_uniform), replace = false)
    bg = [sampl_bg; uni_bg]
    return bg
end

function buffer_and_rasterize(x; to, buffer = 2)
    buffered_points = GO.buffer.(x, buffer)
    rasterize(last, buffered_points; to, fill = true, missingval = false)
end

struct _ContinuousBoyceIndex 
    n_bins::Integer
    bin_overlap::AbstractFloat
    min::Union{AbstractFloat, Nothing}
    max::Union{AbstractFloat, Nothing}
    cor::Function
end

cbi(; n_bins = 101, bin_overlap = 0.1, min = nothing, max = nothing, cor = StatsBase.corspearman) = _ContinuousBoyceIndex(n_bins, bin_overlap, min, max, cor)

## CBI
using StatisticalMeasures
import CategoricalDistributions: classes, pdf
function (m::_ContinuousBoyceIndex)(ŷ::AbstractArray, y)
    positive_class = classes(first(ŷ))|> last
    scores = pdf.(ŷ, positive_class)
    ma = isnothing(m.max) ? maximum(scores) : m.max
    mi = isnothing(m.min) ? minimum(scores) : m.min
    binwidth = m.bin_overlap * (ma - mi)

    return _cbi(scores, y, positive_class, m.n_bins, binwidth, ma, mi, m.cor)

end

const ContinuousBoyceIndex(args...) = _ContinuousBoyceIndex(args...) |> robust_measure |> fussy_measure


function _cbi(scores, y, positive_class, nbins, binwidth, ma, mi, cor)
    binstarts = range(mi, stop=ma-binwidth, length=nbins)
    binends = range(mi + binwidth, stop=ma, length=nbins)

    sorted_indices = sortperm(scores)
    sorted_scores = view(scores, sorted_indices)
    sorted_y = view(y, sorted_indices)

    tot_positive = count(==(positive_class), y)
    tot_negative = length(y) - tot_positive

    n_positive = zeros(Int, nbins)
    n_negative = zeros(Int, nbins)

    @inbounds for i in 1:nbins
        bin_index_first = searchsortedfirst(sorted_scores, binstarts[i])
        bin_index_last = searchsortedlast(sorted_scores, binends[i])
        @inbounds for j in bin_index_first:bin_index_last
            if sorted_y[j] == positive_class 
                n_positive[i] += 1
            end
        end
        n_negative[i] = bin_index_last - bin_index_first + 1 - n_positive[i]
    end

    n_total = n_positive .+ n_negative

    # omit bins with no negative - we don't want to divide by zero
    no_obs = n_negative .== 0
    deleteat!(n_positive, no_obs)
    deleteat!(n_negative, no_obs)
    binstarts = binstarts[.!no_obs]

    binmeans = (n_positive ./ tot_positive) ./ (n_negative ./ tot_negative)
    r = cor(binmeans, binstarts)
    if isnan(r)
        @show binmeans
        @show binstarts
        @show no_obs
        error()
    end
    return r
end

StatisticalMeasures.@trait(_ContinuousBoyceIndex,
       consumes_multiple_observations=true,
       observation_scitype = Finite{2},
       kind_of_proxy=StatisticalMeasures.LearnAPI.Distribution(),
       orientation=Score(),
       external_aggregation_mode=Mean(),
       human_name = "continuous boyce index",
)





#=
import Meshes
function spatial_resample(geometries, y; nfolds, sides = 5, n_tries = 50, rng = Xoshiro(0))
    all_indices = 1:length(geometries)
    ps = Meshes.PointSet(Tuple.(geometries))
    parts = Meshes.partition(ps, Meshes.BlockPartition(sides))
    nparts_presence = sum(p -> any(view(y, p.inds)), parts)
    nparts_absence = sum(p -> !all(view(y, p.inds)), parts)

    nparts_presence > nfolds || error("Trying to divide into $(nfolds) groups, but only $nparts_presence parts have presences")
    nparts_absence > nfolds || error("Trying to divide into $(nfolds) groups, but only $nparts_absence parts have absences")

    for i in 1:n_tries
        parts_groups = rand(rng, 1:nfolds, length(parts))
        train_test_rows = map(1:nfolds) do g
            test_indices = Int64[]
            for i in 1:length(parts)
                if parts_groups[i] == g
                    test_indices = [test_indices; parts[i].inds]
                end
            end
            (Base.setdiff(all_indices, test_indices), test_indices)
        end

        resample_ok = map(train_test_rows) do (train, test)
            Base.sort(Base.unique(y[train])) == [false, true] && 
            Base.sort(Base.unique(y[test])) == [false, true] 
        end |> all

        if resample_ok
            return train_test_rows
        end
    end
    error("Could not find a valid resample after $n_tries tries.")
end
=#