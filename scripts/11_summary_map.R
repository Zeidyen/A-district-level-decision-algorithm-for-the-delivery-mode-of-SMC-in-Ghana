# Summary map: districts by how robust their classification is across the 12 step 4
# variants (data/derived/district_robustness.csv, from the scenario runs).
suppressMessages({library(sf)}); sf_use_s2(FALSE)
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE); adm1$region <- sub(" Region$", "", adm1$shapeName)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- adm1$region[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)]
adm2 <- adm2[adm2$region %in% NORTH, ]
rb <- read.csv("data/derived/district_robustness.csv")
grp <- ifelse(rb$red == rb$n, "Red in all variants",
       ifelse(rb$red >= 9, sprintf("Red in at least 9 of %d", rb$n[1]),
       ifelse(rb$green == 0, "Never green", "Could switch under some assumptions")))
adm2$grp <- grp[match(adm2$shapeName, rb$district)]
lev <- c("Could switch under some assumptions", "Never green", sprintf("Red in at least 9 of %d", rb$n[1]), "Red in all variants")
cols <- c("#7fbf7b", "#f6d58e", "#e08a6e", "#a61c1c")
n_var <- rb$n[1]
nr <- adm1[adm1$region %in% NORTH, ]
png("out/summary_district_map.png", width = 1500, height = 1650, res = 200)
par(mar = c(0.5, 0.5, 3, 0.5))
plot(st_geometry(adm2), col = cols[match(adm2$grp, lev)], border = "white", lwd = 0.6,
     main = sprintf("Can routine delivery replace the SMC campaign?\nDistrict class across %d model variants", n_var),
     cex.main = 0.95)
plot(st_geometry(nr), add = TRUE, lwd = 1.6)
lab <- st_coordinates(st_point_on_surface(nr))
text(lab[, 1], lab[, 2], toupper(nr$region), cex = 0.6, font = 2, col = "grey10")
legend("bottomleft", sprintf("%s (%d)", lev, as.integer(table(factor(adm2$grp, lev)))),
       fill = cols, border = NA, bty = "n", cex = 0.72)
dev.off()
