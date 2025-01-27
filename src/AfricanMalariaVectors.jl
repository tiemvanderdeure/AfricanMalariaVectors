__precompile__(false)
module AfricanMalariaVectors
    using CSV, Tables, DataFrames, StatsBase, Statistics, RCall, GLM
    using Rasters, RasterDataSources, ArchGDAL, NaturalEarth, LibGEOS, Dates
    import MLJBase, Random

    import GeoInterface as GI
    import GeometryOps as GO
    import CategoricalArrays as CA


    import RasterDataSources.URIs: URI

    export SPECIES, LIMITS, DATES, SSPS
    export GAMClassifier
    export load_occurrences, generate_bg
    export load_bioclim, load_kg_classification, load_worldpop
    export map_axis, as_label, rasters_to_cmap, bivariate_cmap_legend, bivariate_colormap, letter_label
    export Bshap


    const DD = Rasters.DimensionalData

    const SPECIES = [
        :arabiensis
        :coluzzii
        :coustani_sl
        :funestus
        :gambiae
        :hancocki
        :marshalli
        :moucheti_sl
        :nili_sl
        :paludis
        :pharoensis
        :rufipes
        :squamosus
        :wellcomei
    ]

    include("gam.jl")
    include("utils.jl")
    include("vectordata.jl")
    include("roi.jl")
    include("rasterdata.jl")
    include("plots.jl")
    include("shapley.jl")
    include("cbi.jl")
    include("gridcv.jl")
end
