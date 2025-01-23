using RasterDataSources, Rasters, ArchGDAL, Dates
using Rasters: DD
const GCMS = [GFDL_ESM4, IPSL_CM6A_LR, MPI_ESM1_2_HR, MRI_ESM2_0, UKESM1_0_LL]
const SSPS = [SSP126, SSP370]
const DATES = [Date(2055), Date(2085)]

# tbd in RastersRasterDataSources extension
function my_rasstack(T::Type{<:WorldClim{<:Future{BioClim}}}; res, date, kw...)
    ras = Raster(RasterDataSources.getraster(T; res, date); kw...)
    @views Rasters.RasterStack(NamedTuple(Symbol("bio$i") => ras[Band = i] for i in 1:19))
end
cropmask(r, m) = mask(read(crop(r; to = m)); with = m)

function load_bioclim(predictors; aggregate= true, aggregation_factor = 5)
    bio = RasterStack(CHELSA{BioClim}, predictors; lazy = true)
    current = read(crop(bio; to = EXTENT))

    bioclim_future = cat(
        (cat(
            (cat(
                (
                    rebuild(
                            RasterStack(
                                CHELSA{Future{BioClim, CMIP6, GCM, SSP}}, predictors; 
                                date, lazy = true, missingval = nothing
                            );
                            dims = dims(bio)
                        )
                 for GCM in GCMS)...;
                dims = Dim{:gcm}(GCMS)
            ) for SSP in SSPS)...;
            dims = Dim{:ssp}(SSPS)
        ) for date in DATES)...;
        dims = Dim{:date}(DATES)
    )
    future = read(crop(bioclim_future; to = EXTENT))

    if aggregate
        current = Rasters.aggregate(mean, current, aggregation_factor; skipmissingval = true)
        future = Rasters.aggregate(mean, future, (X(aggregation_factor), Y(aggregation_factor)); skipmissingval = true)
    end

    current = map(layers(bio), layers(current)) do b, l
        md = Rasters.metadata(b)
        Float32.(l .* md["scale"] .+ md["offset"])
    end |> RasterStack
    future = map(layers(bio), layers(future)) do b, l
        md = Rasters.metadata(b)
        Float32.(l .* md["scale"] .+ md["offset"])
    end |> RasterStack

    return (; current, future)
end

function load_f_lulc(ssp, date; kw...) 
    ssp = "$(string(SSP126)[1:4])_RCP$(string(SSP126)[5:6])"
    yr = year(date)
    path = joinpath(
        ENV["RASTERDATASOURCES_PATH"], 
        "Global 7-land-types LULC projection dataset under SSPs-RCPs", 
        ssp,
        "global_$(ssp)_$yr.tif"
    )

    Raster(path, name = :lulc; kw...)
end

# From the readme. 7 is permanent snow and ice which there isn't any of
const LULC_TYPES = (
    2 => "Forest", 
    3 => "Grassland", 
    4 => "Barren", 
    5 => "Cropland", 
    6 => "Urban", 
)
function load_lulc(roi)
    # from https://zenodo.org/records/4584775
    path = joinpath(
        ENV["RASTERDATASOURCES_PATH"], 
        "Global 7-land-types LULC projection dataset under SSPs-RCPs",
        "global_LULC_2015.tif"
    )
    isdir(path) || download_lulc() 
    lulc = Raster(path, name = :lulc)
    lulc_res = CA.categorical(resample(lulc; to = roi, method = :mode))
    current = CA.recode(lulc_res, missing, LULC_TYPES...) # default to missing so permanent water is missing

    lulc_f = cat(
            (cat(
            (
                resample(load_f_lulc(ssp, date); to = roi)
            for ssp in SSPS)...;
            dims = Dim{:ssp}(SSPS)
        ) for date in DATES)...;
        dims = Dim{:date}(DATES)
    )
    future = CA.recode(CA.categorical(lulc_f), missing, LULC_TYPES...)

    return (; current, future)
end

function load_worldpop(to)
    dir = joinpath(RasterDataSources.rasterpath(), "WorldPop")
    if !isdir(dir) mkdir(dir) end
    filepath = joinpath(dir, "ppp_2020_1km_Aggregated.tif")
    fileurl = URI("https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/0_Mosaicked/ppp_2020_1km_Aggregated.tif")
    maybe_download(fileurl, filepath)
    pop = Raster(filepath; lazy = true) |> replace_missing
    return resample(read(crop(pop; to)); to, method = :sum) # aggregate might be better here?
end

function download_lulc() 
    url = URI("https://zenodo.org/records/4584775/files/Global%207-land-types%20LULC%20projection%20dataset%20under%20SSPs-RCPs.zip")
    dir = joinpath(RasterDataSources.rasterpath(), "Global 7-land-types LULC projection dataset under SSPs-RCPs.zip")
    maybe_download(url, dir)
    destpath = joinpath(RasterDataSources.rasterpath(), "Global 7-land-types LULC projection dataset under SSPs-RCPs")
    isdir(destpath) || mkdir(destpath)
    try 
        run(`unzip $dir -d $destpath`)
    catch
        error("Could not unzip file at $dir, please unzip manually")
    end

end

## Faster broadcasting over categorical arrays
function Base.Broadcast.broadcasted(::typeof(==), ca::CA.CategoricalArray, x)
    index = CA.DataAPI.invrefpool(ca)
    if haskey(index, x)
        CA.refs(ca) .== index[x]
    else
        broadcast(x -> ismissing(x) ? missing : false, ca)
    end
end
Base.Broadcast.broadcasted(::typeof(==), ca::CA.CategoricalArray, x::Missing) = fill(missing, size(ca))
