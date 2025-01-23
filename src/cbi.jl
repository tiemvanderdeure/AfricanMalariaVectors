## Define the Continuous Boyce Index
using StatisticalMeasures
import CategoricalDistributions: classes, pdf
struct _ContinuousBoyceIndex 
    n_bins::Integer
    bin_overlap::AbstractFloat
    min::Union{AbstractFloat, Nothing}
    max::Union{AbstractFloat, Nothing}
    cor::Function
end

cbi(; n_bins = 101, bin_overlap = 0.1, min = nothing, max = nothing, cor = StatsBase.corspearman) = _ContinuousBoyceIndex(n_bins, bin_overlap, min, max, cor)

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
    isnan(r) && error()
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
