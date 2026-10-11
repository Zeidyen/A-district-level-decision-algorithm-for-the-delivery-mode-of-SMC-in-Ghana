# Covariates for the stacked machine-learning model of routine-contact slips (25_ml_stack.py).
#
# All layers on the 1 km Ghana grid of the step 1 travel-time surfaces (tt_walk_any_gha.tif), so
# the learners can be trained on all Ghana DHS clusters and predict every northern grid cell.
#
#   tt_walk, tt_motor   travel time to the nearest public facility (step 1; log2, floor 5 min)
#   logpopdens          population density (log)
#   fac_5km, chps_5km   public facilities / CHPS within about 5 km (9 x 9 cell window; log1p)
#   viirs               night-time lights, VIIRS 2016 (WorldPop covariate; log1p of the 1 km mean)
#   dst_road, dst_roadx distance to the nearest major road / road intersection, OSM 2016 (log1p km)
#   dst_water           distance to the nearest waterway, OSM 2016 (log1p km)
#   dst_urban           distance to artificial surfaces, ESA CCI 2015 (class 190; signed log km)
#   dst_crop            distance to rainfed cropland, ESA CCI 2015 (class 011; signed log km)
#   dst_tree            distance to tree cover, ESA CCI 2015 (class 040; signed log km)
#   slope, elevation    SRTM slope (WorldPop) and elevation
#   rain_annual         mean annual rainfall 2020-22, CHIRPS (mm)
#   rain_wetshare       share of annual rain falling in June-September
# WorldPop 100 m layers are averaged to 1 km (missing distances = inside the feature = 0). JRC surface water covers only the north, so it is
# not used here.
#
# Cluster values are under-5-weighted means over the DHS displacement buffer (2 km urban, 5 km
# rural), as in 03_dhs_clusters.R. Northern clusters also get their district (as 22_district_cv.R).
#
# Outputs: data/derived/ml_covariates_gha.tif, ml_cluster_covariates.csv, ml_grid_covariates.csv

source("scripts/lib_inla.R")
tmpl <- rast("data/derived/tt_walk_any_gha.tif")
wp <- function(f) {
  r <- rast(file.path("data/raw/worldpop_cov", f))
  resample(r, tmpl, method = "average")
}
# WorldPop distance layers are in km. The ESA CCI ones are signed (negative inside the land-cover
# class), so all distances get a signed log, sign(x) log(1 + |x|); a missing distance is 0. Edge
# gaps in the other layers (a few coastal or border cells) are filled from the 5 x 5 neighbourhood.
km <- function(r) { r <- ifel(is.na(r), 0, r); sign(r) * log1p(abs(r)) }
fill <- function(r) cover(r, focal(r, w = 5, fun = "mean", na.rm = TRUE))

ttw <- rast("data/derived/tt_walk_any_gha.tif"); ttm <- rast("data/derived/tt_motor_any_gha.tif")
fac <- read.csv("data/derived/facilities.csv")
cnt <- function(sel) {
  r <- rasterize(cbind(fac$lon[sel], fac$lat[sel]), tmpl, fun = "count", background = 0)
  log1p(focal(r, w = matrix(1, 9, 9), fun = "sum", na.rm = TRUE))
}
rain <- rast("data/derived/chirps_monthly_gha.tif")
rain <- rain[[grep("^202[0-2]-", names(rain))]]
mon <- as.integer(substr(names(rain), 6, 7))
rain_annual <- resample(sum(rain) / 3, tmpl, method = "bilinear")
rain_wet <- resample(sum(rain[[mon %in% 6:9]]) / sum(rain), tmpl, method = "bilinear")
elev <- resample(rast("data/raw/elevation/elevation/GHA_elv_msk.tif"), tmpl, method = "bilinear")

S <- c(log2(max(ttw, 5)), log2(max(ttm, 5)), rast("data/derived/logpopdens_gha.tif"),
       cnt(rep(TRUE, nrow(fac))), cnt(fac$is_chps %in% c(TRUE, "TRUE")),
       fill(log1p(wp("gha_viirs_100m_2016.tif"))),
       km(wp("gha_osm_dst_road_100m_2016.tif")), km(wp("gha_osm_dst_roadintersec_100m_2016.tif")),
       km(wp("gha_osm_dst_waterway_100m_2016.tif")),
       km(wp("gha_esaccilc_dst190_100m_2015.tif")), km(wp("gha_esaccilc_dst011_100m_2015.tif")),
       km(wp("gha_esaccilc_dst040_100m_2015.tif")),
       fill(wp("gha_srtm_slope_100m.tif")), fill(elev), rain_annual, rain_wet)
names(S) <- c("tt_walk", "tt_motor", "logpopdens", "fac_5km", "chps_5km", "viirs", "dst_road", "dst_roadx",
              "dst_water", "dst_urban", "dst_crop", "dst_tree", "slope", "elevation", "rain_annual", "rain_wetshare")
S <- mask(S, tmpl)
writeRaster(S, "data/derived/ml_covariates_gha.tif", overwrite = TRUE)

# ---- clusters: under-5-weighted buffer means -----------------------------------------
ge <- st_read("data/raw/dhs/GHGE8AFL/GHGE8AFL.shp", quiet = TRUE)
ge <- ge[ge$LATNUM != 0, ]
buf <- vect(st_buffer(ge, ifelse(ge$URBAN_RURA == "U", 2000, 5000) / 111320))
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") + rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
u5 <- resample(u5, tmpl); u5[is.na(u5) | u5 < 0] <- 0
num <- extract(S * u5, buf, fun = sum, na.rm = TRUE)[, -1]
den <- extract(u5, buf, fun = sum, na.rm = TRUE)[, 2]
plain <- extract(S, buf, fun = mean, na.rm = TRUE)[, -1]
M <- as.matrix(num) / den
M[den == 0, ] <- as.matrix(plain)[den == 0, ]
cv <- as.data.frame(M)
names(cv) <- names(S)
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)])
adm2 <- adm2[adm2$region %in% NORTH, ]
dist <- adm2$shapeName[sapply(st_intersects(ge, adm2), function(i) if (length(i)) i[1] else NA)]
out <- data.frame(cluster = ge$DHSCLUST, district = dist, cv)
write.csv(out, "data/derived/ml_cluster_covariates.csv", row.names = FALSE)
cat("clusters:", nrow(out), "| any missing covariate:", sum(!complete.cases(out[, names(S)])), "\n")

# ---- northern grid cells (same cells as north_grid(), used by step 4) ------------------
gd <- north_grid()
gv <- extract(S, gd[, c("x", "y")])[, -1]
write.csv(data.frame(cell = gd$cell, x = gd$x, y = gd$y, gv), "data/derived/ml_grid_covariates.csv", row.names = FALSE)
cat("grid cells:", nrow(gd), "| with any missing covariate:", sum(!complete.cases(gv)), "\n")
print(round(cor(out[, names(S)], use = "pairwise"), 2)[1:6, ])
