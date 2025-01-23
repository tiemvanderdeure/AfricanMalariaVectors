function thin(x, cutoff; geometrycolumn = :geometry)
    if !(x isa AbstractVector{<:GI.NamedTuplePoint}) && Tables.istable(x)
        geoms = Tables.getcolumn(Tables.Columns(x), geometrycolumn)
        indices = _thin(geoms, cutoff)
        return Tables.subset(x, indices)
    elseif x isa AbstractVector
        _geoms = unique(geoms)
        indices = _thin(_geoms, cutoff)
        return _geoms[indices]
    else   
        ArgumentError("Supply a Table with a geometry column, or an iterable of points")
    end
end

function _thin(geoms::AbstractVector, cutoff)
    for i in 1:length(geoms)
        GI.geomtrait(geoms[i]) == GI.PointTrait() && ArgumentError("$(geoms[i]) is not a point")
    end

    dist_matrix = pairwise(Haversine(), geoms)

    dist_mask = dist_matrix .< cutoff
    dist_sum = vec(sum(dist_mask; dims = 1))
    s = size(dist_sum, 1)

    indices = collect(1:s)

    for i in 1:s
        m = maximum(dist_sum)
        if m == 1
            break
        else
            drop = rand(findall(dist_sum .== m))
            dist_sum .-= @view dist_mask[drop, :]
            drop_indices = 1:s .!= drop
            dist_mask = @view dist_mask[drop_indices, drop_indices]
            dist_sum = @view dist_sum[drop_indices]
            s -= 1
            indices = @view indices[drop_indices]
        end
    end
    return indices
end
