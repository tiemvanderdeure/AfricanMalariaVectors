using Shapley, Combinatorics, Tables, Statistics, ProgressMeter
import SpeciesDistributionModels as SDM

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

#=
n = 200

x = (a = rand(n), b = rand(n), c= rand(n), d = rand(n), e = rand(n), f = rand(n))
Z = (a = x.a .+ 0.5, b = x.b .+ 0.1, c = x.c, d = x.d, e = x.e, f = rand(n))
pred(x) = (x.a .+ x.b .+ x.c).^x.d

shapvals = Bshap(pred, x, Z)

pred(Z) .- pred(X)

shapvals.e |> sum
### scales as 2^N (length(collect(combinations(1:N))))

combs = 
length(collect(combs))
l = length(comb)
nfeatures - 1 - l

for comb in vcat([[]], collect(combs))
    l = length(comb)
    n_permutations = factorial(l) * factorial(nfeatures - 1 - l)
    for c in comb
        X_[c] .= Z[c]
    end
    pred_base = pred(X_)
end

@time merge(X_[(:a, :b)], Z[(:c,)])

@time for c in comb
    X_[c] .= Z[c]
end

collect(combs)

# particular case of Shapley
perms = permutations(1:2)

changes = NamedTuple(k => zeros(length(perms), Tables.rowcount(X)) for k in keys(tab))
tab = deepcopy(X)

map(enumerate(perms)) do (i, p)
    for k in keys(tab)
        tab[k] .= X[k]
    end
    refvalue = pred(tab)
    for j in p
        tab[j] .= Z[j]
        newvalue = pred(tab)
        changes[j][i,:] .= (newvalue .- refvalue)
        refvalue = newvalue
    end
end

map(x -> mean(x; dims = 2), changes)

X = (a = rand(2), b = rand(2))
Z = (a = X.a .+ 0.5, b = X.b .- 0.1)

shapley(x -> x.a .+ x.b, Shapley.MonteCarlo(10_000), X, 1, Z)
=#