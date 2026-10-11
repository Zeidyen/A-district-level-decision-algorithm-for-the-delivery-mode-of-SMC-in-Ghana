# Step 1: travel time from every ~1 km pixel to the nearest routine platform,
# across the five northern SMC regions.
#
# Facilities: data/derived/facilities.csv from 01b_facilities.R. Since 10 Oct 2026 the
# Ghana Health Service public facility finder (1,658 CHPS in the five northern regions);
# before that, the Maina et al. 2019 list (~212 northern CHPS).
#
# Facilities up to 50 km outside the study area are kept, so children near the
# boundary can reach a facility across it. Facilities in Burkina Faso and Togo
# are not included (cross-border care is ignored).
#
# Outputs (data/derived/): tt_{walk,motor}_{any,chps}.tif in minutes,
# district_access_u5.csv (under-5 children by travel-time band and district),
# and out/access_map.png.

suppressMessages({
  library(terra); library(sf); library(gdistance); library(raster); library(readxl)
})
sf_use_s2(FALSE)

STUDY_REGIONS <- c("Northern Region", "Savannah Region", "North East Region",
                   "Upper East Region", "Upper West Region")
BANDS <- c(0, 30, 60, 120, Inf)        # minutes
BAND_LABELS <- c("<30 min", "30-60 min", "60-120 min", ">120 min")

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
study <- adm1[adm1$shapeName %in% STUDY_REGIONS, ]
stopifnot(nrow(study) == 5)

# Districts take the region that contains most of their area.
inside <- st_intersects(st_point_on_surface(adm2), adm1)
adm2$region <- adm1$shapeName[sapply(inside, `[`, 1)]
districts <- adm2[adm2$region %in% STUDY_REGIONS, ]

# Facilities: the cleaned file from 01b_facilities.R (GHS facility finder by default; set
# FACILITY_SOURCE=maina there for the earlier Maina et al. list).
fac <- read.csv("data/derived/facilities.csv")
cat("facility source:", fac$source[1], "\n")
fac <- st_as_sf(fac, coords = c("lon", "lat"), crs = 4326)
buffer <- st_buffer(st_union(study), 0.45)   # ~50 km in degrees
fac <- fac[st_intersects(fac, buffer, sparse = FALSE)[, 1], ]
cat(nrow(fac), "facilities in study area + buffer,", sum(fac$is_chps), "of them CHPS\n")

# Travel time on the MAP friction surfaces (minutes per metre).
extent_box <- ext(vect(buffer))
travel_time <- function(friction_file, points) {
  fr <- crop(rast(friction_file), extent_box)
  tr <- transition(raster(1 / fr), function(x) mean(x), directions = 8)
  tr <- geoCorrection(tr, type = "c")
  rast(accCost(tr, st_coordinates(points)))
}

dir.create("data/derived", showWarnings = FALSE)
tt <- list()
for (mode in c("walk", "motor")) {
  ff <- sprintf("data/raw/map/friction_%s_2019.tif",
                c(walk = "walking", motor = "motorised")[mode])
  for (target in c("any", "chps")) {
    pts <- if (target == "any") fac else fac[fac$is_chps, ]
    r <- travel_time(ff, pts)
    r <- mask(crop(r, vect(study)), vect(study))
    names(r) <- paste0("tt_", mode, "_", target)
    writeRaster(r, sprintf("data/derived/tt_%s_%s.tif", mode, target), overwrite = TRUE)
    tt[[names(r)]] <- r
  }
}

# Sanity check against MAP's published travel time to healthcare (which uses a
# larger, different facility set, so ours should be similar or longer).
for (mode in c("walk", "motor")) {
  mp <- rast(sprintf("data/raw/map/map_tt_%s_healthcare.tif",
                     c(walk = "walking", motor = "motorised")[mode]))
  mp <- mask(crop(mp, vect(study)), vect(study))
  ours <- resample(tt[[paste0("tt_", mode, "_any")]], mp)
  v <- na.omit(cbind(values(ours), values(mp)))
  cat(sprintf("%s: median ours %.0f min vs MAP %.0f min; correlation %.2f\n",
              mode, median(v[, 1]), median(v[, 2]), cor(v[, 1], v[, 2], method = "spearman")))
}

# Under-5 children by travel-time band and district.
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") +
      rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
u5 <- crop(u5, vect(study))
u5[is.na(u5) | u5 < 0] <- 0

dv <- vect(districts)
rows <- list()
for (nm in names(tt)) {
  t <- resample(tt[[nm]], u5, method = "bilinear")
  band <- classify(t, cbind(BANDS[-length(BANDS)], BANDS[-1], seq_along(BAND_LABELS)),
                   right = FALSE)
  for (b in seq_along(BAND_LABELS)) {
    pop <- extract(u5 * (band == b), dv, fun = sum, na.rm = TRUE)[, 2]
    rows[[length(rows) + 1]] <- data.frame(
      region = districts$region, district = districts$shapeName,
      measure = nm, band = BAND_LABELS[b], u5 = round(pop))
  }
}
out <- do.call(rbind, rows)
write.csv(out, "data/derived/district_access_u5.csv", row.names = FALSE)

# Regional summary printed for the worklog.
reg <- aggregate(u5 ~ region + measure + band, out, sum)
reg$share <- ave(reg$u5, reg$region, reg$measure, FUN = function(x) round(100 * x / sum(x)))
reg$band <- factor(reg$band, BAND_LABELS)
print(reg[order(reg$measure, reg$region, reg$band), ], row.names = FALSE)

# Map: walking time to any public facility, and to listed CHPS.
dir.create("out", showWarnings = FALSE)
png("out/access_map.png", width = 2000, height = 1100, res = 150)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 5))
pal <- hcl.colors(4, "YlOrRd")
fac_in <- st_intersects(fac, st_union(study), sparse = FALSE)[, 1]
for (nm in c("tt_walk_any", "tt_walk_chps")) {
  plot(tt[[nm]], breaks = BANDS, col = pal, axes = FALSE, box = FALSE,
       main = c(tt_walk_any = "Walking time to nearest public facility",
                tt_walk_chps = "Walking time to nearest listed CHPS")[nm],
       plg = list(legend = BAND_LABELS, title = "minutes"))
  plot(st_geometry(districts), add = TRUE, border = "grey40", lwd = 0.3)
  plot(st_geometry(study), add = TRUE, lwd = 1.2)
  # Buffer facilities shape travel times near the border but are not drawn.
  shown <- fac_in & (if (nm == "tt_walk_chps") fac$is_chps else TRUE)
  plot(st_geometry(fac[shown, ]), add = TRUE, pch = 16, cex = 0.25)
}
dev.off()
cat("done\n")
