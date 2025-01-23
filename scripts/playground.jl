using Rasters, CSV, DataFrames, ArchGDAL
using CairoMakie
malariadata = CSV.read("../malariadata/00 Africa 1900-2015 SSA PR database (260617).csv", DataFrame)
filter!(r -> r.AREA_TYPE == "Point" && r.YY > 1980 && r.YY <= 2010, malariadata)

p_gambiae = Raster("preds_24sep/current/gambiae/raster.tif") |> replace_missing

extr = extract(p_gambiae, malariadata, geometrycolumn = (:Long, :Lat), skipmissing = false)
malariadata.p_gambiae .= getindex.(extr, 2)
malaria_gambiae = dropmissing(malariadata, [:p_gambiae, Symbol("PfPR2-10")])
corspearman(malaria_gambiae[!, Symbol("PfPR2-10")], malaria_gambiae.p_gambiae)


scatter(
    malaria_gambiae[!, Symbol("PfPR2-10")], malaria_gambiae.p_gambiae,
    color = PfPR2-10)







#### playing around starts here
ensembles_tmax = map(occs_thinned, bgs) do x, bg
    data = SDM.sdmdata(x, bg; predictors)
    SDM.sdm(data, models)
end
ensembles_tmin = map(occs_thinned, bgs) do x, bg
    data = SDM.sdmdata(x, bg; predictors = (:bio1, :bio6, :bio7, :bio12, :bio15, :bio14))
    SDM.sdm(data, models)
end

bg = reduce(hcat, map(b -> collect(b[PREDICTORS]), AnophelesSDMs.sample_rstack(bioclim.current, 1000; wts = bm)))
bg = reduce(hcat, map(b -> collect(b[PREDICTORS]), all_occs))

cor(bg')

predictors

pr_gamb = map(b -> SDM.predict(ensembles.gambiae, b; reducer = mean), bioclim)
pr_gamb_tmax = map(b -> SDM.predict(ensembles_tmax.gambiae, b; reducer = mean), bioclim)
pr_gamb_tmin = map(b -> SDM.predict(ensembles_tmin.gambiae, b; reducer = mean), bioclim)

using CairoMakie
fig = Figure()
for (i, p) in enumerate(pr_gamb)
    ax = map_axis(fig[1, i])
    plot!(dropdims(pr_gamb[i]; dims = (3,4,5)); colorrange = (0,1))
    ax = map_axis(fig[2, i])
    plot!(dropdims(pr_gamb_tmax[i]; dims = (3,4,5)); colorrange = (0,1))
    ax = map_axis(fig[3, i])
    plot!(dropdims(pr_gamb_tmin[i]; dims = (3,4,5)); colorrange = (0,1))
end
fig

plot(dropdims(pr_gamb_tmin.future .- pr_gamb.future; dims = (3,4,5)); colorrange = (-0.2, 0.2), colormap = :vik)


pr_gamb.future |> skipmissing |> mean

pr_gamb[2]

### playing around ends here