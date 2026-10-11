# Environmental layers for the five northern SMC regions, for inspection.
#   maps:  mean annual rainfall 2020-22 (CHIRPS), share of annual rain in Jun-Sep,
#          seasonal surface water share (JRC GSW 2021), elevation (SRTM 30 arc-s),
#          under-5 confirmed malaria incidence, annual and Jun-Oct (DHIMS2 2020-24, by district)
#   chart: monthly rainfall by region, 2020-22 mean, with the SMC cycles marked
# Outputs: out/env_maps.png, out/env_rain_months.png

suppressMessages({library(terra); library(sf); library(geodata)})
sf_use_s2(FALSE)
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE); adm1$region <- sub(" Region$", "", adm1$shapeName)
nr <- adm1[adm1$region %in% NORTH, ]; v <- vect(nr)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
clip <- function(r) mask(crop(r, v), v)

rain <- rast("data/derived/chirps_monthly_gha.tif")
rain <- rain[[grep("^202[0-2]-", names(rain))]]
mon <- as.integer(substr(names(rain), 6, 7))
annual <- clip(sum(rain) / 3)
wetshare <- clip(sum(rain[[mon %in% 6:9]]) / sum(rain))
water <- clip(rast("data/derived/jrc_seasonal_share.tif"))
dem_file <- "data/raw/elevation/GHA_elv_msk.tif"
elev <- clip(if (file.exists(dem_file)) rast(dem_file) else elevation_30s("GHA", path = "data/raw/elevation"))
bd <- read.csv("data/derived/district_burden.csv")         # 09_burden.R (under-5, DHIMS)
adm2n <- adm2; adm2n$region <- adm1$region[sapply(st_intersects(st_point_on_surface(adm2n), adm1), `[`, 1)]
adm2n <- adm2n[adm2n$region %in% NORTH, ]
adm2n$annual <- bd$annual_inc_u5[match(adm2n$shapeName, bd$district)]
adm2n$season <- bd$season_inc_u5[match(adm2n$shapeName, bd$district)]
choro <- function(x, title, unit) {
  br <- pretty(x, 6); pal <- hcl.colors(length(br) - 1, "Reds", rev = TRUE)
  plot(st_geometry(adm2n), col = pal[cut(x, br, include.lowest = TRUE)], border = "white", lwd = 0.4,
       main = title, cex.main = 0.95)
  plot(st_geometry(nr), add = TRUE, lwd = 0.9)
  legend("bottomleft", levels(cut(x, br, include.lowest = TRUE)), fill = pal, border = NA, bty = "n", cex = 0.7, title = unit)
}

panel <- function(r, title, pal, unit, rng = NULL) {
  plot(r, col = pal, axes = FALSE, box = FALSE, main = title, cex.main = 0.95,
       range = rng, plg = list(title = unit, cex = 0.8))
  plot(st_geometry(adm2), add = TRUE, border = "grey55", lwd = 0.2)
  plot(st_geometry(nr), add = TRUE, lwd = 0.9)
}
png("out/env_maps.png", width = 2400, height = 1650, res = 170)
par(mfrow = c(2, 3), mar = c(0.5, 0.5, 2.5, 4.5))
panel(annual, "Mean annual rainfall, 2020-22 (CHIRPS)", hcl.colors(50, "YlGnBu", rev = TRUE), "mm")
panel(wetshare, "Share of annual rain falling Jun-Sep (SMC season)", hcl.colors(50, "Purples", rev = TRUE), "share")
panel(clamp(water, 0, 0.5), "Seasonal surface water share (JRC 2021)", hcl.colors(50, "Blues", rev = TRUE),
      "share (capped 0.5)")
panel(elev, "Elevation (SRTM)", hcl.colors(50, "Terrain"), "m")
choro(adm2n$annual, "Under-5 confirmed malaria per child per year (DHIMS 2020-24)", "cases / child / yr")
choro(adm2n$season, "Under-5 confirmed malaria per child, Jun-Oct (DHIMS 2020-24)", "cases / child")
dev.off()

# Monthly rainfall by region (area mean, 2020-22 average).
rm <- extract(rain, v, fun = mean, na.rm = TRUE)[, -1]
clim <- t(sapply(seq_len(nrow(nr)), function(i) tapply(unlist(rm[i, ]), mon, mean)))
rownames(clim) <- nr$region
clim <- clim[NORTH, ]
cols <- c("#1b6ca8", "#3a9e5f", "#c77c11", "#8e44ad", "#c0392b")
png("out/env_rain_months.png", width = 1700, height = 950, res = 180)
par(mar = c(4, 4.5, 3, 9), xpd = FALSE)
plot(NA, xlim = c(1, 12), ylim = c(0, max(clim) * 1.08), xaxt = "n", las = 1,
     xlab = "", ylab = "Mean monthly rainfall, 2020-22 (mm)",
     main = "Rain peaks in Jul-Sep in every northern region; SMC cycles run Jun-Sep", cex.main = 0.95)
rect(5.5, 0, 9.5, max(clim) * 1.08, col = "grey92", border = NA)
text(7.5, max(clim) * 1.04, "SMC cycles (Jun-Sep)", cex = 0.8, col = "grey30")
axis(1, 1:12, month.abb)
for (i in seq_len(nrow(clim))) lines(1:12, clim[i, ], col = cols[i], lwd = 2.2, type = "o", pch = 16, cex = 0.6)
par(xpd = TRUE)
legend("topright", inset = c(-0.27, 0), rownames(clim), col = cols, lwd = 2.2, bty = "n", cex = 0.85)
dev.off()
print(round(clim))
cat("done\n")
