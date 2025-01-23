using RCall, DataFrames, GLM
import MLJModelInterface as MMI
import CategoricalArrays as CA
R"library(mgcv)"

mutable struct GAMClassifier <:  MMI.Probabilistic
    k::Int
    method::Symbol
    gamma::Float64
end
GAMClassifier(; k, method = :REML, gamma = 1.0) = GAMClassifier(k, method, gamma)
function MMI.fit(m::GAMClassifier, verbosity, X, y)
    df = DataFrame(X)
    schema = MMI.schema(X)
    cont_vars = schema.names[findall(.!(schema.scitypes .<: MMI.Finite))]
    cat_vars = schema.names[findall(schema.scitypes .<: MMI.Finite)]
    @assert !in(:y, schema.names)
    y_boolean = Bool.(MMI.int(y) .- 1)
    df.y = y_boolean
    decode = MMI.classes(y)

    formjoin = if !isempty(cont_vars) && !isempty(cat_vars) " - 1 + " else "" end

    formula = "y~" * join(["s($x, k = $(m.k), bs = 'cr')" for x in cont_vars], "+") * formjoin * join(["as.factor($x)" for x in cat_vars], "+")

    gam = R"mgcv::gam(as.formula($formula), family = 'binomial', data = $df, drop.unused.levels = FALSE, method = $(string(m.method)), gamma = $(m.gamma))"
    # Get all the necessary object from the R object here so that predict can be threaded
    coefs = rcopy(R"coefficients($gam)")
    coefnames = rcopy(R"names(coefficients($gam))")

    smooths = map(eachindex(cont_vars)) do i     
        smooth = R"$gam$smooth[[$i]]"
        xk = rcopy(R"$smooth$xp")
        F = rcopy(R"$smooth$F")
        Q = rcopy(R"qr.Q(attr($smooth, 'qrc'), complete = TRUE)")
        (xk, F, Q)
    end
    gam = (;gam, coefs, coefnames, smooths)

    cache = nothing
    report = nothing
    fitresult = (gam, (cont_vars, cat_vars), decode)
    return (fitresult, cache, report)
end

function MMI.predict(m::GAMClassifier, (gam, vars, decode), Xnew)
    cont_vars, cat_vars = vars
    nrow = Tables.rowcount(Xnew)

    p = zeros(nrow)
    M = zeros(nrow, m.k)

    for (i, v) in enumerate(cont_vars)
        xk, F, Q = gam.smooths[i]
        Fm = reshape(F, m.k, m.k)
        X = Tables.getcolumn(Xnew, v)
        cubregspline!.(eachrow(M), X; xk, Fm)
        QC = Q* [0; gam.coefs[startswith("s($v)").(gam.coefnames)]]
        p .+= M * QC
    end

    for (i, v) in enumerate(cat_vars)
        X = Tables.getcolumn(Xnew, v)
        c = gam.coefs[[findfirst(==("as.factor($v)$l"), gam.coefnames) for l in levels(X)]]
        p .+= view(c, CA.refs(X))
    end

    intercept_idx = findfirst(==("(Intercept)"), gam.coefnames)
    !isnothing(intercept_idx) && (p .+= gam.coefs[intercept_idx])

    p .= GLM.linkinv.(LogitLink(), p)
    MMI.UnivariateFinite(decode, [1 .- p, p])
end

function Rpredict(mach, Xnew)
    df = DataFrame(Xnew)
    gam = mach.fitresult[1].gam
    p = rcopy(R"predict($gam, $df)")
    p .= GLM.linkinv.(LogitLink(), p)
    MMI.UnivariateFinite(decode, [1 .- p, p])   
end

### Adapted from the MGCV source code
function cubregspline!(Xp, x; xk::AbstractVector, Fm::Matrix{Float64})
    nk = length(xk)
    kmin = first(xk)
    kmax = last(xk)
    if x <= kmin
        j = 1
        h = xk[2] - kmin
        xik = x - kmin
        cjm = -xik * h / 3
        cjp = -xik * h / 6
        ajm = 1 - xik / h
        ajp = xik / h
    elseif x >= kmax
        j = nk-1
        h = kmax - xk[end-1]
        xik = x - kmax
        cjm = xik * h / 6
        cjp = xik * h / 3
        ajm = - xik / h
        ajp = 1 + xik / h
    else
        j = searchsortedfirst(xk, x)  - 1
        h = xk[j+1] - xk[j]
        ajm = (xk[j+1] - x)
        ajp = (x-xk[j])
        cjm = (ajm*(ajm*ajm/h - h))/6;
        cjp = (ajp*(ajp*ajp/h - h))/6;
        ajm /= h
        ajp /= h
    end
    @inbounds @views Xp .= cjm .* Fm[:, j] .+ cjp .* Fm[:, j+1]
    Xp[j] += ajm
    Xp[j+1] += ajp
    return
end
cubregspline(x; xk::AbstractVector, Fm::Matrix{Float64}) = cubregspline!(similar(xk), x; xk, Fm)

