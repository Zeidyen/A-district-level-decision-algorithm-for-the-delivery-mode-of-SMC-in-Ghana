# Ghana-wide travel time to the nearest public facility (walking and motorised),
# so that clusters from every region can inform the step 2 covariate effects.
# Same method as 02_access.R (which stays the source of the northern district
# tables). Facilities: data/derived/facilities.csv (01b_facilities.R), all of Ghana.
#
# Also writes log population density on the same grid, the second prediction
# covariate (a raster is needed for every covariate the map predicts from, so
# cluster-only variables such as wealth cannot be used in the step 2 model).
#
# Outputs (data/derived/): tt_{walk,motor}_any_gha.tif, logpopdens_gha.tif

suppressMessages({
  library(terra); library(sf); library(gdistance); library(raster); library(readxl)
})
sf_use_s2(FALSE)

gha <- st_union(st_read("data/raw/gha_ADM1.geojson", quiet = TRUE))
fac <- read.csv("data/derived/facilities.csv")             # from 01b_facilities.R
fac <- st_as_sf(fac, coords = c("lon", "lat"), crs = 4326)
cat(nrow(fac), "public facilities in Ghana\n")

box <- ext(vect(st_buffer(gha, 0.1)))
for (mode in c("walk", "motor")) {
  fr <- crop(rast(sprintf("data/raw/map/friction_%s_2019.tif",
                          c(walk = "walking", motor = "motorised")[mode])), box)
  tr <- geoCorrection(transition(raster(1 / fr), function(x) mean(x), 8), type = "c")
  r <- mask(rast(accCost(tr, st_coordinates(fac))), vect(gha))
  names(r) <- paste0("tt_", mode, "_any")
  writeRaster(r, sprintf("data/derived/tt_%s_any_gha.tif", mode), overwrite = TRUE)
  cat(mode, ": median", round(global(r, median, na.rm = TRUE)[[1]]), "min\n")
}

# Population density (all ages, people per km2) from the WorldPop grid, put on
# the travel-time grid. Under-5 counts are a fixed share of all ages at this
# scale, so they serve as the density proxy.
tt <- rast("data/derived/tt_walk_any_gha.tif")
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") +
      rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
u5[is.na(u5) | u5 < 0] <- 0
dens <- resample(u5, tt, method = "bilinear") / cellSize(tt, unit = "km")
lpd <- mask(log1p(dens), tt)
names(lpd) <- "logpopdens"
writeRaster(lpd, "data/derived/logpopdens_gha.tif", overwrite = TRUE)
cat("done\n")
