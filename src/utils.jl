## Background point
function generate_bg(rng, x, all_occs, bioc, bm = Rasters.boolmask(bioc); fraction_uniform, fraction_sampling)
    n_p = length(x)
    other_occs = Base.setdiff(all_occs, x)
    uni_bg = Rasters.sample(rng, bioc, floor(Int, n_p * fraction_sampling); skipmissing = true, geometry = (X,Y))
    sampl_bg = sample(rng, other_occs, floor(Int, n_p * fraction_uniform), replace = false)
    bg = [sampl_bg; uni_bg]
    return bg
end

## Throw errors if a file is not found
function check_file(path, url)
    isfile(path) || error(
        "
            Expected file at $path, but it does not exist.
            This file is available at $url
            Please download it manually and store it at the specified path
        "
    )
    return
end

## Downloads a file if it doesn't exist already
function maybedownload(url, path; force = false)
    if force || !isfile(path)
        @info "Starting download for $url"
        download(url, path)
    end
end

using RCall
## Installs R packages
function Rrequire(pkg)
    only(rcopy(R"""require($pkg, character.only = T)""")) && return nothing
    @info "$pkg not available in R, installing now...."
    reval("""install.packages("$pkg", method="wget", repos="https://cloud.r-project.org")""")
    reval("""library("$pkg", character.only = T)""")
  
    return nothing
end
