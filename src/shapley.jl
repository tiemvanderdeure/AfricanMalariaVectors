using Shapley, Combinatorics, Tables, Statistics, ProgressMeter
import SpeciesDistributionModels as SDM

# adapt from shapley library to handle rasterstack inputs
function shapley(alg::Shapley.Algorithm, e::SDM.SDMensemble, X::AbstractRasterStack)
    bm = boolmask(X)
    preds = Tables.columntable(X[bm])
    preds = map(x -> nonmissingtype(eltype(x)).(x), preds)
    outras = rebuild(X; data = 
        map(_ -> Array{Union{Missing, Float64}}(missing, size(X)), layers(X))
    )
    
    f(x) = SDM.predict(e, x; reducer = mean)

    values = Shapley.shapley(f, alg, preds)
    
    for k in keys(outras)
        outras[k][bm] .= values[k]
    end
    return outras
end

# implement Bshap as defined by http://proceedings.mlr.press/v119/sundararajan20b/sundararajan20b.pdf
function Bshap(predict, X::AbstractRasterStack, Z::AbstractRasterStack; kw...)
    @assert dims(X) == dims(Z)
    bm = boolmask(X) .&& boolmask(Z)
    outras = rebuild(X; data = 
        map(_ -> Array{Union{Missing, Float64}}(missing, size(X)), layers(X))
    )
    values = Bshap(predict, map.(identity, X[bm]), map.(identity, Z[bm]); kw...)
    for k in keys(outras)
        outras[k][bm] .= values[k]
    end
    return outras
end
function Bshap(predict, X, Z; features = nothing)
    X_ = Tables.columntable(X)
    Z_ = Tables.columntable(Z)

    nfeatures = length(X_)
    features = isnothing(features) ? Tables.columnnames(X_) : Tables.columnnames(X_)[features]
    shapvalues = NamedTuple(i => zeros(Tables.rowcount(X_)) for i in features)
    cache = Dict{Set{Symbol}, Vector{Float64}}()

    for i in 0:(nfeatures-1)
        combs = Combinatorics.combinations(features, i)
        for comb in combs
            l = length(comb)
            # number of permutations per f starting with comb
            n_permutations = factorial(l) * factorial(nfeatures - 1 - l)
            data = merge(X_, Z_[Tuple(comb)])
            set = Set(comb)
            if haskey(cache, set) 
                baseline = cache[set]
            else # this should only be the case for i == 0
                baseline = predict(data)
                cache[set] = baseline
            end
            Threads.@threads for f in features
                if !(f ∈ comb)
                    newset = Set(vcat(comb, f))
                    if haskey(cache, newset)
                        new = cache[newset]
                    else
                        newdata = merge(data, Z_[(f,)])
                        new = predict(newdata)
                        cache[newset] = new
                    end
                    shapvalues[f] .+= (new .- baseline) .* n_permutations
                end
            end
            # delete the set from the cache to save memory
            delete!(cache, set)
        end
    end

    for k in keys(shapvalues)
        shapvalues[k] ./= factorial(nfeatures)
    end
    return shapvalues
end
