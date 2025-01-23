species_metadata() = CSV.read("data/species.csv", DataFrame; types = Dict(:species_level => Bool))

function load_occurrences()
    # Malaria Atlas Project data
    mapdata = CSV.read("data/map_data.csv", DataFrame; missingstring = "NA")
    map_species = CSV.read("data/map_species.csv", DataFrame; delim = ',') |> dropmissing!
    dropmissing!(mapdata, [:longitude, :latitude, :year_end]) # allow missing years for this dataset??
    mapdata.geometry = map((x,y) -> (x,y), mapdata.longitude, mapdata.latitude)
    rename!(mapdata, :species => :map_species, :citation => :source,)
    mapdata = innerjoin(mapdata, map_species; on = :map_species)
    mapdata = mapdata[!, [:year_start, :year_end, :species, :geometry]]

    # Massey et al 2016 data
    massey = CSV.read("data/Bionomics Africa.csv", DataFrame)
    massey_species = CSV.read("data/massey_species.csv", DataFrame) |> dropmissing!
    dropmissing!(massey, [:long, :lat, :year_end])
    filter!(x -> x.area_type == "point", massey)
    massey.geometry = map((x,y) -> (x,y), massey.long, massey.lat)
    massey.source = massey.citation
    rename!(massey, :species => :massey_species)
    massey = innerjoin(massey, massey_species; on = :massey_species)
    massey = massey[!, [:year_start, :year_end, :species, :geometry]]

    # Snow / Kyalo et al 2017 dataset
    vectors = CSV.read("data/dataverse_files/Africa Vectors database_1898-2016.csv", DataFrame)
    snow_species = CSV.read("data/snow_species.csv", DataFrame) |> dropmissing!
    dropmissing!(vectors, [:Long, :Lat])
    vectors.geometry = map((x,y) -> (x,y), vectors.Long, vectors.Lat)
    species_cols = startswith.(names(vectors), "An") .|| startswith.(names(vectors), "SS")
    vectors_melt = DataFrames.stack(vectors, species_cols, [:geometry, :YeEnd, :YeStart], variable_name = :snow_species)
    dropmissing!(vectors_melt)

    other_sib_species = mapreduce(vcat, filter(row -> !ismissing(row["Other sib species names"]), eachrow(vectors))) do row
        species = strip.(split(row["Other sib species names"], (',', ';')))
        DataFrame(; snow_species = species, row.geometry, row.YeEnd, row.YeStart, value = "Y")
    end

    vectors_melt = [vectors_melt; other_sib_species]
    rename!(vectors_melt, "YeStart" => :year_start, "YeEnd" => :year_end)
    snow = innerjoin(vectors_melt, snow_species; on = :snow_species)
    DataFrames.select!(snow, Not([:value, :snow_species]))

    all_data = [mapdata; massey; snow]
    unique!(all_data)
    filter!(x -> x.year_start < 2010 && x.year_end > 1980, all_data)

    sp_md = species_metadata()
    all_data = innerjoin(all_data, sp_md; on = :species)
    return all_data
end

function load_malaria(; dir = "data")
    path = joinpath(dir, "00 Africa 1900-2015 SSA PR database (260617).csv")
    isfile(path) || error(
        "
            Expected malaria database to be at $path, but it does not exist.
            The malaria database is available at https://doi.org/10.7910/DVN/Z29FR0
            Please download it manually and store it at the specified path
        "
    )
    CSV.read(path, DataFrame)
end