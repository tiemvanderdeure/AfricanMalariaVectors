## Shared dimensions
const GCMS = Dim{:gcm}([GFDL_ESM4, IPSL_CM6A_LR, MPI_ESM1_2_HR, MRI_ESM2_0, UKESM1_0_LL])
const SSPS = Dim{:ssp}([SSP126, SSP370])
const DATES = Ti([Date(2055), Date(2085)])

## Climate data

function load_bioclim(predictors; aggregate= true, aggregation_factor = 5)
    current = load_current_bioclim(predictors, aggregate, aggregation_factor)

    future = (@d load_future_bioclim.(Ref(predictors), GCMS, SSPS, DATES, aggregate, aggregation_factor)) |>
        RasterSeries |> Rasters.combine

    return (; current, future)
end

function load_current_bioclim(predictors, aggregate, aggregation_factor)
    bio = RasterStack(CHELSA{BioClimPlus}, predictors; lazy = true, missingval=nothing)
    current = read(crop(bio; to = EXTENT))
    # Fix an issue with CHELSA data where 0s are typemax
    map((:gst, :gsp), (1e3, 1e6)) do k, v
        if haskey(current, k)
            current[k][current[k] .> v] .= 0 # CHELSA data is stored as UInt32, so we need to add one to the values
        end
    end

    if aggregate
        current = Rasters.aggregate(mean, current, aggregation_factor; skipmissingval = true)
    end
    return current
end

function load_future_bioclim(predictors, gcm, ssp, date, aggregate, aggregation_factor)
    data = RasterStack(
        CHELSA{Future{BioClimPlus, CMIP6, gcm, ssp}}, predictors; 
        date, lazy = true, missingval = nothing
    )
    data = read(crop(data; to = EXTENT))
    if aggregate 
        data = Rasters.aggregate(mean, data, aggregation_factor; skipmissingval = true)
    end
    return data
end

### Land use data

# LULC types from the readme of the publication. 
# 7 is permanent snow and ice which there isn't any of
# 1 is water wich we recode to missing
const LULC_TYPES = (
    2 => "Forest", 
    3 => "Grassland", 
    4 => "Barren", 
    5 => "Cropland", 
    6 => "Urban", 
)

function load_lulc(roi)
    # from https://zenodo.org/records/4584775
    basepath = maybe_download_lulc() 
    lulc = load_current_lulc(basepath)
    lulc_res = CA.categorical(resample(lulc; to = roi, method = :mode))
    current = CA.recode(lulc_res, missing, LULC_TYPES...) # default to missing so permanent water is missing

    lulc_f = (@d load_f_lulc.(basepath, SSPS, DATES)) |> RasterSeries |> Rasters.combine
    lulc_f_res = CA.categorical(resample(lulc_f; to = roi, method = :mode))
    future = CA.recode(lulc_f_res, missing, LULC_TYPES...)

    return (; current, future)
end

load_current_lulc(basepath = maybe_download_lulc()) = Raster(joinpath(basepath, "global_LULC_2015.tif"), name = :lulc)

function load_f_lulc(basepath, ssp, date; kw...) 
    ssp = "$(string(ssp)[1:4])_RCP$(string(ssp)[5:6])"
    yr = year(date)
    path = joinpath(basepath, ssp, "global_$(ssp)_$yr.tif")
    Raster(path, name = :lulc; kw...)
end

function maybe_download_lulc() 
    destpath = joinpath(RasterDataSources.rasterpath(), "Global 7-land-types LULC projection dataset under SSPs-RCPs")
    # Download and un-zip
    if !isdir(destpath)
        mkdir(destpath)
        url = URI("https://zenodo.org/records/4584775/files/Global%207-land-types%20LULC%20projection%20dataset%20under%20SSPs-RCPs.zip")
        dir = joinpath(RasterDataSources.rasterpath(), "Global 7-land-types LULC projection dataset under SSPs-RCPs.zip")
        maybe_download(url, dir)
        try 
            run(`unzip $dir -d $destpath`)
        catch
            error("Could not unzip file at $dir, please unzip manually")
        end
    end
    return destpath
end

### Population data

function load_worldpop(roi)
    current = load_current_worldpop(roi)
    future = load_future_worldpop(current)
    return (; current, future)
end

function load_current_worldpop(to)
    dir = joinpath(RasterDataSources.rasterpath(), "WorldPop")
    if !isdir(dir) mkdir(dir) end
    filepath = joinpath(dir, "ppp_2020_1km_Aggregated.tif")
    fileurl = URI("https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/0_Mosaicked/ppp_2020_1km_Aggregated.tif")
    maybe_download(fileurl, filepath)
    pop = Raster(filepath; lazy = true) |> replace_missing
    return resample(read(crop(pop; to)); to, method = :sum)
end

function get_pop_ssp(ssp)
    dir = joinpath(RasterDataSources.rasterpath(), "WIC_Population_Projections")
    if !isdir(dir) mkdir(dir) end
    ssp = string(ssp)[1:4]
    filename = "PROJresult_AGE_$(ssp)_V14.csv"
    url = URI("https://zenodo.org/records/14718294/files/$filename?download=1")
    maybe_download(url, joinpath("data", filename))
    return CSV.read(joinpath("data", filename), DataFrame)
end


function load_future_worldpop(current)
    countries = naturalearth("ne_10m_admin_0_countries_iso") |> DataFrame
    countries.ISO_N3 = parse.(Int, countries.ISO_N3)
    yrs = Tuple(year.(DATES))

    pop_changes = map(SSPS) do ssp
        pop_ssp = get_pop_ssp(ssp)

        filter!(x -> x.agest >= 0, pop_ssp)
        pop_future_ssp1 = filter(x -> x.Time in (2020, yrs...), pop_ssp)
        pop_totals_ssp1 = DataFrames.combine(groupby(pop_future_ssp1, [:region, :Time]), :pop => sum => :pop)
        pop_unstack = unstack(pop_totals_ssp1, :Time, :pop)
        for yr in Symbol.(yrs)
            pop_unstack[!, yr] = pop_unstack[!, yr] ./ pop_unstack[!, Symbol(2020)]
        end
        pop_unstack.ISO_N3 = parse.(Int, replace.(pop_unstack.region, "reg" => ""))
        pop_with_geometries = innerjoin(pop_unstack, countries, on = :ISO_N3)
        pop_change_rasterized = rasterize(last, pop_with_geometries; to = current, fill = Symbol.(yrs))
        cat(layers(pop_change_rasterized)...; dims = DATES)
    end
    relative_pop = cat(pop_changes...; dims = SSPS)
    relative_pop = replace_missing(relative_pop, 1) # such that pixels with missing data do not completely disappear
    return relative_pop .* current
end
