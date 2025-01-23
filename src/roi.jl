const LIMITS = Float64[-18, 52, -35, 25]
const EXTENT = Rasters.Extents.Extent(X = (LIMITS[1], LIMITS[2]), Y = (LIMITS[3], LIMITS[4]))

function get_land()
    countries = naturalearth("ne_10m_admin_0_countries")
    intersecting_countries = filter(g -> GO.coverage(g, LIMITS...) > 0, countries.geometry)
    countries_multipoly = reduce((x,y) -> GO.union(GO.GEOS(), x, y; target = GO.MultiPolygonTrait()), intersecting_countries);
    filter(g -> GO.coverage(g, LIMITS...) > 0, collect(GI.getpolygon(countries_multipoly))) |> GI.MultiPolygon
end

function get_africa()
    countries = naturalearth("ne_10m_admin_0_countries")
    africa_geoms = countries.geometry[countries.CONTINENT .== "Africa"]
    reduce(LibGEOS.union, africa_geoms)
end

function get_roi(land = get_land(), bio12 = load_bioclim((:bio12,)).current, africa = get_africa())
    africa_land = LibGEOS.intersection(land, africa)
    notdesert = replace_missing(bio12 .>= 200, false) # desert has <200 mm precipitation
    notdesert_multipoly = mypolygonize(notdesert) # polygonize
    main_poly = argmax(GO.area, GI.getpolygon(notdesert_multipoly)) # find the main area to exclude tiny spots
    roi = GO.buffer(main_poly, 7) # buffer to include +- half of the Sahara
    roi_ras = rasterize(roi; to = notdesert, missingval = false, fill = true)
    roi_ras[Y = -Inf .. 20] .= true # include all the rest of Africa (madagascar etc)
    mask(roi_ras; with = africa_land, boundary = :inside) # mask out sea and lakes
end

function mypolygonize(r::Raster)
    xs, ys = map(lookup(r)) do l
        range(first(l), last(l); length = length(l))
    end

    GO.polygonize(xs, ys, parent(r)) # polygonize
end

