# Crop the CHIRPS v2.0 monthly Africa rasters (Nov 2019 - Jan 2023, the span of
# the DHS card dates) to Ghana and stack them. Raw files: data/raw/chirps/.
suppressMessages(library(terra))
f <- sort(list.files("data/raw/chirps", pattern = "tif$", full.names = TRUE))
g <- vect("data/raw/gha_ADM1.geojson")
r <- crop(rast(f), ext(g) + 0.2)
r[r < 0] <- NA
names(r) <- sub(".*v2\\.0\\.([0-9]{4})\\.([0-9]{2})\\.tif$", "\\1-\\2", f)
writeRaster(r, "data/derived/chirps_monthly_gha.tif", overwrite = TRUE)
ue <- g[g$shapeName == "Upper East Region"]
cat("Upper East mean rainfall (mm), 2020:\n")
print(round(unlist(extract(r[[grep("^2020", names(r))]], ue, fun = mean, na.rm = TRUE)[, -1])))
