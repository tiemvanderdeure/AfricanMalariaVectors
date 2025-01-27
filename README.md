# Climate change favours African malaria vectors

This repository contains code for "Climate change favours African malaria vectors" by van der Deure, Nogués-Bravo, Njotto & Stensgaard, which is currently in preparation. In this publication, we investigate the current and future climate and land use conditions for six dominant and widely spread Africa malaria vectors.

All code used in the analysis is contained in this repository, and all data is either automatically downloaded when the scripts are run, or can be downloaded using the instructions below (the code will prompt when the file is not found in the expected location).

## Understanding this repository
This repository is formatted as a Julia package. The `Project.toml` specifies which packages it relies on. The `src` folder contains its internal code and the `data` folder contain data files. The only code you have to interact with directly to reproduce all analysis in the paper is located in the `scripts` folder. 

## Reproducing the analysis
First clone this repository and navigate to the folder.

### Installing julia and instantiating the environment
All analysis is performed in the Julia 1.11. The recommended way to install Julia is through [juliaup](https://github.com/JuliaLang/juliaup). 

After installing julia, navigate to this repository, then type `]` in the Julia REPL to enter package mode, and type `activate .` followed by `instantiate` build a reproducible environment that contains all packages necessary to run the code.

This project uses R to fit GAMs and to download Malaria Atlas Project data, interfacing between the languages using RCall.jl. A working R installation is required and will be automatically detected and the necessary R packages installed. R version 4.4.2 was used to run the code.

### Downloading data
Most of the data will be automatically downloaded. This includes large amounts (10s of GB) of raster data, which will be downloaded to a path specified in `ENV["RASTERDATASOURCES_PATH]`. You can set this to any path (e.g. a hard-drive) by running `ENV["RASTERDATASOURCES_PATH"] = "my/data/path"` in the Julia REPL. The data will then be stored in several sub-folders. Vector and malaria observation data is read from the `data` folder.

Two sources cannot be directly downloaded and need to be downloaded manually:
- The geo-coded inventory of anophelines by Snow (2017). It is available at https://doi.org/10.7910/DVN/NQ6CUN and contains a .csv file that should be copied to `data`, so it can be loaded from `data\Africa Vectors database_1898-2016.csv`.
- The Plasmodium falciparum prevalence datast by Snow (2017). It is available at https://doi.org/10.7910/DVN/Z29FR0 and contains a .csv file that should be copied to `data`, so it can be loaded from `data\00 Africa 1900-2015 SSA PR database (260617).csv`.

### Running the code
The `scripts` folder contains three Julia scripts: `scripts/main.jl` runs the analysis, while `scripts/figures.jl` and `scripts/tables.jl` use the Julia objects to create figures and tables in the publications. These are saved in the `images` and `tables` folders. 

You could run these from the command line by navigating the this folder, opening a Julia REPL by running the `julia` command, and then running each of the three the scripts using an include statements (`include("scripts/main.jl")` etc).

Alternatively, you can use an environment like Visual Studio Code to run them line-by-line.

This code makes use of multi-threading to speed up computation, see the [relevant section of the Julia manual](https://docs.julialang.org/en/v1/manual/multi-threading/) for how to run Julia with multiple threads.

## Contact
For any questions about the code in this repository, contact Tiem van der Deure at email tvd@sund.ku.dk.

Making an issue in this repository works as well.
