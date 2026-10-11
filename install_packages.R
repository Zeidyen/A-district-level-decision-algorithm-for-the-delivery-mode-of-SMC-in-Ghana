# R packages used by the pipeline (developed with R 4.5.1). INLA comes from its own repository.
pkgs <- c("dplyr", "tidyr", "sf", "terra", "raster", "gdistance", "geodata", "malariaAtlas", "haven",
          "readxl", "glmmTMB", "MASS")
install.packages(setdiff(pkgs, rownames(installed.packages())))
if (!requireNamespace("INLA", quietly = TRUE))
  install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
