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
    map(1:resampler.nfolds) do f
        test = rows[fold .== f]
        train = Base.setdiff(rows, test)
        return (train, test)
    end
end

