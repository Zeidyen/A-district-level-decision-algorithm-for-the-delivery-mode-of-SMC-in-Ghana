# Terrain (SRTM 30 arc-s elevation and slope) in the north, and whether it adds
# anything to routine-contact timeliness beyond travel time.
# Cluster values are under-5-weighted means over the DHS displacement buffer (as 03).
# Outputs: out/terrain_maps.png; prints the cluster-level check.
suppressMessages({library(terra); library(sf); library(dplyr); library(glmmTMB)})
sf_use_s2(FALSE)
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE); adm1$region <- sub(" Region$", "", adm1$shapeName)
nr <- adm1[adm1$region %in% NORTH, ]; v <- vect(nr)
elev <- rast(list.files("data/raw/elevation", "GHA_elv_msk.tif$", recursive = TRUE, full.names = TRUE)[1])
slope <- terrain(elev, "slope", unit = "degrees")

png("out/terrain_maps.png", width = 2200, height = 1000, res = 170)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 4.5))
plot(mask(crop(elev, v), v), col = hcl.colors(50, "Terrain"), axes = FALSE, box = FALSE,
     main = "Elevation (m, SRTM ~1 km)", plg = list(title = "m"))
plot(st_geometry(nr), add = TRUE, lwd = 0.9)
plot(clamp(mask(crop(slope, v), v), 0, 6), col = hcl.colors(50, "YlOrBr", rev = TRUE), axes = FALSE, box = FALSE,
     main = "Slope (degrees, capped at 6)", plg = list(title = "deg"))
plot(st_geometry(nr), add = TRUE, lwd = 0.9)
dev.off()

sv <- values(mask(crop(slope, v), v)); sv <- sv[!is.na(sv)]
cat(sprintf("North slope: median %.2f deg, 95th pct %.2f, %% cells > 5 deg: %.1f\n",
            median(sv), quantile(sv, .95), 100 * mean(sv > 5)))

# Cluster check: does terrain predict per-contact slip beyond motorised travel time?
cl <- read.csv("data/derived/dhs_clusters.csv") |>
  inner_join(read.csv("data/derived/dhs_clusters_timeliness.csv"), by = "cluster") |>
  filter(contacts > 0, region %in% NORTH)
pts <- st_as_sf(cl, coords = c("LONGNUM", "LATNUM"), crs = 4326)
buf <- vect(st_buffer(pts, ifelse(cl$urban == "U", 2000, 5000) / 111320))
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") + rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
u5 <- resample(u5, elev); u5[is.na(u5)] <- 0
wm <- function(r) { e <- extract(c(r * u5, u5), buf, fun = sum, na.rm = TRUE); e[, 2] / e[, 3] }
cl$elev <- wm(elev); cl$slope <- wm(slope)
cl$q <- cl$slip56 / cl$contacts
sp <- function(a, b) round(cor(a, b, method = "spearman", use = "complete"), 2)
cat(sprintf("Northern clusters %d | Spearman with q: elevation %.2f, slope %.2f | slope vs motor time %.2f\n",
            nrow(cl), sp(cl$elev, cl$q), sp(cl$slope, cl$q), sp(cl$slope, cl$tt_motor_any)))
cl$lt <- log2(pmax(cl$tt_motor_any, 5)); cl$id <- seq_len(nrow(cl))
or <- function(m, t) { b <- fixef(m)$cond[t]; s <- sqrt(vcov(m)$cond[t, t])
  sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * s), exp(b + 1.96 * s)) }
m0 <- glmmTMB(cbind(slip56, contacts - slip56) ~ lt + (1 | id), data = cl, family = binomial)
m1 <- glmmTMB(cbind(slip56, contacts - slip56) ~ lt + scale(slope) + scale(elev) + (1 | id), data = cl, family = binomial)
cat("OR per doubling of motor time:", or(m0, "lt"), "-> with terrain", or(m1, "lt"), "\n")
cat("OR per SD slope:", or(m1, "scale(slope)"), "| per SD elevation:", or(m1, "scale(elev)"),
    "| AIC", round(AIC(m0), 1), "->", round(AIC(m1), 1), "\n")

# Is elevation a regional proxy? Refit with region fixed effects.
m2 <- glmmTMB(cbind(slip56, contacts - slip56) ~ lt + scale(slope) + scale(elev) + region + (1 | id),
              data = cl, family = binomial)
cat("With region effects -> OR per SD elevation:", or(m2, "scale(elev)"), "| slope:", or(m2, "scale(slope)"),
    "| motor time:", or(m2, "lt"), "| AIC", round(AIC(m2), 1), "\n")
cat("Mean cluster elevation by region (m):",
    paste(names(tapply(cl$elev, cl$region, mean)), round(tapply(cl$elev, cl$region, mean)), collapse = " | "), "\n")
