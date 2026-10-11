# Rainy-season (SMC season, Jun-Sep) travel time to the nearest public facility.
#
# MAP's friction surfaces are season-free. Here cells that flood seasonally are
# slowed in the wet season, using JRC Global Surface Water v1.4 (30 m; Pekel et al.
# 2016, Nature) "seasonality" for 2021: months with water present, 1-11 = seasonal,
# 12 = permanent. The 30 m layer is averaged to the 1 km MAP grid as the share of
# each cell under seasonal water.
#
# Rule: in the wet season a cell whose seasonal-water share >= F takes the friction
# MAP itself gives open water (0.06 min/m walking, ~1 km/h; checked: the median
# walking friction of cells >= 50% permanent water). Road cells (MAP motorised
# friction below walking friction) stay at their dry value, i.e. a bridge or
# causeway is assumed, unless roads are set to flood too.
#   main:   F = 0.10, roads passable
#   mild:   F = 0.25, roads passable
#   severe: F = 0.05, roads flood as well (unpaved tracks cut)
#   barrier: as severe, but flooded cells near-impassable (1 min/m, ~17 h per km):
#            an extreme bound, not a realistic scenario
# Facilities as in 02_access.R (Maina list, north + 50 km, (0,0) rows dropped).
#
# Outputs: data/derived/tt_{walk,motor}_wet_{main,mild,severe}.tif (north + buffer),
#          data/derived/wet_access_summary.csv, out/wet_access_map.png

suppressMessages({library(terra); library(sf); library(gdistance); library(raster); library(readxl)})
sf_use_s2(FALSE)

NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
WATER_FRICTION <- 0.06
RULES <- list(main = list(F = 0.10, roads_flood = FALSE),
              mild = list(F = 0.25, roads_flood = FALSE),
              severe = list(F = 0.05, roads_flood = TRUE),
              barrier = list(F = 0.05, roads_flood = TRUE, friction = 1))   # near-impassable bound

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE); adm1$region <- sub(" Region$", "", adm1$shapeName)
study <- adm1[adm1$region %in% NORTH, ]
buffer <- st_buffer(st_union(study), 0.45)
E <- ext(vect(buffer))

fac <- read_excel("data/raw/facilities/ssa_public_facilities.xlsx")
fac <- fac[fac$Country == "Ghana" & !(fac$Lat == 0 & fac$Long == 0), ]
fac <- st_as_sf(fac, coords = c("Long", "Lat"), crs = 4326)
fac <- fac[st_intersects(fac, buffer, sparse = FALSE)[, 1], ]

fw <- crop(rast("data/raw/map/friction_walking_2019.tif"), E)
fm <- crop(rast("data/raw/map/friction_motorised_2019.tif"), E)
seas_file <- "data/derived/jrc_seasonal_share.tif"
if (!file.exists(seas_file)) {
  s <- crop(vrt(list.files("data/raw/jrc_gsw", "^seasonality", full.names = TRUE)), E)
  writeRaster(resample(s >= 1 & s <= 11, fw, method = "average"), seas_file, overwrite = TRUE)
}
seas <- rast(seas_file); seas[is.na(seas)] <- 0
road <- fm < fw

travel_time <- function(fr, pts) {
  tr <- geoCorrection(transition(raster(1 / fr), function(x) mean(x), directions = 8), type = "c")
  rast(accCost(tr, st_coordinates(pts)))
}

tt <- list(dry_walk = travel_time(fw, fac), dry_motor = travel_time(fm, fac))
for (nm in names(RULES)) {
  r <- RULES[[nm]]
  flood <- seas >= r$F & (if (r$roads_flood) TRUE else !road)
  wf <- if (is.null(r$friction)) WATER_FRICTION else r$friction
  fw2 <- ifel(flood, max(fw, wf), fw)                      # never make a cell faster
  fm2 <- ifel(flood, max(fm, wf), fm)
  cat(sprintf("%s: %.1f%% of cells flood\n", nm, 100 * global(flood, "mean", na.rm = TRUE)[[1]]))
  tt[[paste0(nm, "_walk")]] <- travel_time(fw2, fac)
  tt[[paste0(nm, "_motor")]] <- travel_time(fm2, fac)
  for (mode in c("walk", "motor"))
    writeRaster(tt[[paste0(nm, "_", mode)]], sprintf("data/derived/tt_%s_wet_%s.tif", mode, nm), overwrite = TRUE)
}

# ---- children affected -----------------------------------------------------------
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") * 0.75 +
      rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
u5 <- resample(u5, fw); u5[is.na(u5) | u5 < 0] <- 0
reg <- rasterize(vect(study), fw, field = "region")
rv <- as.data.frame(c(reg, u5), na.rm = FALSE); names(rv) <- c("region", "u5")

rows <- list()
for (nm in names(RULES)) for (mode in c("walk", "motor")) {
  d <- values(tt[[paste0("dry_", mode)]])[, 1]; w <- values(tt[[paste0(nm, "_", mode)]])[, 1]
  for (rg in c(NORTH, "All north")) {
    i <- !is.na(rv$region) & !is.na(d) & !is.na(w) & (rg == "All north" | rv$region == rg)
    k <- rv$u5[i]
    rows[[length(rows) + 1]] <- data.frame(rule = nm, mode = mode, region = rg,
      dry_mean_min = round(sum(d[i] * k) / sum(k)), wet_mean_min = round(sum(w[i] * k) / sum(k)),
      pct_children_plus30 = round(100 * sum(k[(w[i] - d[i]) > 30]) / sum(k), 1),
      pct_children_x1.5 = round(100 * sum(k[w[i] > 1.5 * pmax(d[i], 5)]) / sum(k), 1))
  }
}
summ <- do.call(rbind, rows)
write.csv(summ, "data/derived/wet_access_summary.csv", row.names = FALSE)
print(summ[summ$mode == "walk", ], row.names = FALSE)
cat("\nMotorised, all north:\n"); print(summ[summ$mode == "motor" & summ$region == "All north", ], row.names = FALSE)

# ---- map: extra walking minutes in the wet season (main rule) -------------------------
dl <- mask(crop(tt$main_walk - tt$dry_walk, vect(study)), vect(study))
png("out/wet_access_map.png", width = 1600, height = 1500, res = 170)
par(mar = c(1, 1, 3, 5))
plot(clamp(dl, 0, 120), col = hcl.colors(50, "Blues", rev = TRUE), axes = FALSE, box = FALSE,
     main = "Extra walking time to nearest facility in the rainy season (min, main rule)",
     cex.main = 0.9)
plot(st_geometry(st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)), add = TRUE, border = "grey60", lwd = 0.3)
plot(st_geometry(study), add = TRUE, lwd = 1.2)
dev.off()
cat("done\n")
