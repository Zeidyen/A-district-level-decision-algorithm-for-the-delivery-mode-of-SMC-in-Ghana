# Study site map: Ghana's 16 regions, the five SMC regions of the north shaded, their 55 districts
# outlined, and an inset locating Ghana in Africa. Greyscale, for print.
#
# Inputs: data/raw/gha_ADM1.geojson, gha_ADM2.geojson (geoBoundaries gbOpen); world boundaries from
# the geodata package (downloaded to data/raw/world/ on first run). Output: out/study_site_map.png

suppressMessages({ library(sf); library(ggplot2); library(geodata) })
sf_use_s2(FALSE)
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm1$region <- sub(" Region$", "", adm1$shapeName)
adm1$smc <- adm1$region %in% NORTH
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
pt <- suppressWarnings(st_point_on_surface(adm2))
adm2$region <- adm1$region[as.integer(st_intersects(pt, adm1))]
dist_n <- adm2[adm2$region %in% NORTH, ]
cat("SMC districts outlined:", nrow(dist_n), "\n")

lab <- suppressWarnings(st_point_on_surface(adm1))
lab <- cbind(adm1[, c("region", "smc")] |> st_drop_geometry(), st_coordinates(lab))
lab$Y[lab$region == "Upper East"] <- lab$Y[lab$region == "Upper East"] + 0.05
lab[lab$region == "North East", c("X", "Y")] <- c(-0.35, 10.42)   # clear of the Northern boundary
lab$name <- ifelse(lab$smc, toupper(lab$region), lab$region)

world <- world(resolution = 3, level = 0, path = "data/raw/world")
world <- st_as_sf(world)
africa <- st_crop(world, xmin = -20, xmax = 55, ymin = -36, ymax = 38)
gha <- st_union(adm1)

main <- ggplot() +
  geom_sf(data = adm1, aes(fill = smc), colour = "grey35", linewidth = 0.35) +
  geom_sf(data = dist_n, fill = NA, colour = "grey55", linewidth = 0.15) +
  geom_sf(data = adm1[adm1$smc, ], fill = NA, colour = "black", linewidth = 0.7) +
  geom_text(data = lab[lab$smc, ], aes(X, Y, label = name), size = 4.6, fontface = "bold") +
  geom_text(data = lab[!lab$smc, ], aes(X, Y, label = name), size = 3.4, colour = "grey25") +
  scale_fill_manual(values = c(`TRUE` = "grey78", `FALSE` = "white"), guide = "none") +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_line(colour = "grey92"), axis.text = element_text(size = 11))

# north arrow and 100 km scale bar (1 degree of longitude at 5.5 N is about 110.7 km)
bb <- st_bbox(gha); x0 <- bb["xmax"] - 0.55; y0 <- bb["ymin"] + 0.35
km100 <- 100 / (111.32 * cos(5.5 * pi / 180))
main <- main +
  annotate("segment", x = x0, xend = x0, y = y0 + 0.5, yend = y0 + 1.0, arrow = arrow(length = unit(0.3, "cm")), linewidth = 0.8) +
  annotate("text", x = x0, y = y0 + 1.18, label = "N", size = 5.5, fontface = "bold") +
  annotate("rect", xmin = x0 - km100, xmax = x0, ymin = y0, ymax = y0 + 0.08, fill = "black") +
  annotate("text", x = x0 - km100 / 2, y = y0 - 0.15, label = "100 km", size = 4.2)

inset <- ggplot() +
  geom_sf(data = africa, fill = "white", colour = "grey60", linewidth = 0.15) +
  geom_sf(data = gha, fill = "black", colour = "black") +
  coord_sf(xlim = c(-18, 52), ylim = c(-35, 37), expand = FALSE) +
  theme_void() + theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6),
                       panel.background = element_rect(fill = "white", colour = NA))

png("out/study_site_map.png", width = 2700, height = 2500, res = 220)
grid::grid.newpage()
print(main, vp = grid::viewport(x = 0.37, y = 0.5, width = 0.74, height = 1))
print(inset, vp = grid::viewport(x = 0.87, y = 0.76, width = 0.23, height = 0.3))
grid::grid.text("Ghana in Africa", x = 0.87, y = 0.93, gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.rect(x = 0.795, y = 0.52, width = 0.03, height = 0.022, gp = grid::gpar(fill = "grey78", col = "black", lwd = 1.5))
grid::grid.text("SMC region", x = 0.818, y = 0.52, just = "left", gp = grid::gpar(fontsize = 13))
grid::grid.lines(x = c(0.78, 0.81), y = 0.48, gp = grid::gpar(col = "grey55", lwd = 1))
grid::grid.text("District boundary", x = 0.818, y = 0.48, just = "left", gp = grid::gpar(fontsize = 13))
grid::grid.lines(x = c(0.78, 0.81), y = 0.44, gp = grid::gpar(col = "grey35", lwd = 1.2))
grid::grid.text("Region boundary", x = 0.818, y = 0.44, just = "left", gp = grid::gpar(fontsize = 13))
grid::grid.text("5 regions, 55 districts", x = 0.78, y = 0.39, just = "left", gp = grid::gpar(fontsize = 13))
dev.off()
cat("done\n")
