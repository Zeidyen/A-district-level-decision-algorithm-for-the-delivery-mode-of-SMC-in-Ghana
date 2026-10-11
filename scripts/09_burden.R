# Step 3: malaria burden in SMC-eligible children (3-59 months) during the months SMC protects:
# doses are given June-September and each protects for about 28 days, so June-October,
# per ~1 km pixel and per district, across the five northern regions. Under-5 only.
#
#   season incidence (district) = DHIMS2 confirmed malaria cases under 5 in Jun-Oct,
#                                 2020-24 mean per year / DHIMS2 under-5 population
#   season cases (pixel)        = district season incidence x WorldPop 2025 children 3-59 m
#
# Observed under-5 cases replace the earlier MAP all-age incidence x age ratio x season share
# (changed 9 Oct 2026: "we are focused on 0-5"). No season share or age conversion is
# needed: the Jun-Oct under-5 cases are counted directly.
#
# Notes:
# - These are cases under the current SMC campaign; step 4 back-calculates the no-SMC
#   baseline as before. They are confirmed cases only (presumptive excluded; under-5 testing
#   rates are 92-100% of suspected cases), and miss children treated outside public
#   facilities, so they understate true burden. Totals scale with that; district classes
#   in step 4 do not depend on burden.
# - Under 5 in DHIMS is 0-59 m; SMC covers 3-59 m. The rate is applied to WorldPop 3-59 m.
# - Rates are uniform within a district (DHIMS is district-level); pixels differ only in
#   the number of children.
#
# Outputs: data/derived/season_incidence_u5.tif, season_cases_u5.tif, district_burden.csv,
#          out/burden_map.png

suppressMessages({library(terra); library(sf); library(dplyr)})
sf_use_s2(FALSE)
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
norm <- function(v) gsub("[^a-z]", "", tolower(v))

dh <- read.csv("data/raw/dhims/dhims_national_2015_2025.csv", check.names = FALSE) |>
  filter(Region %in% NORTH, year %in% 2020:2024)
dist <- dh |> group_by(Region, District, adm2_key) |>
  summarise(season_cases_u5_yr = sum(conf_u5[m %in% 6:10]) / 5,
            annual_cases_u5_yr = sum(conf_u5) / 5,
            pop_u5 = mean(Pop_5), .groups = "drop") |>
  mutate(season_inc = season_cases_u5_yr / pop_u5, annual_inc = annual_cases_u5_yr / pop_u5,
         season_share_u5 = season_cases_u5_yr / annual_cases_u5_yr)
stopifnot(nrow(dist) == 55)

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE); adm1$region <- sub(" Region$", "", adm1$shapeName)
nr <- adm1[adm1$region %in% NORTH, ]
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- adm1$region[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)]
adm2 <- adm2[adm2$region %in% NORTH, ]
adm2$key <- norm(adm2$shapeName)
adm2$season_inc <- dist$season_inc[match(adm2$key, dist$adm2_key)]
stopifnot(!anyNA(adm2$season_inc))

# ---- surfaces -----------------------------------------------------------------
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") * 0.75 +
      rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
u5 <- mask(crop(u5, vect(nr)), vect(nr))
u5[is.na(u5) | u5 < 0] <- 0
season_inc <- rasterize(vect(adm2), u5, field = "season_inc")
season_inc[is.na(season_inc)] <- 0                       # slivers outside any district
cases <- season_inc * u5
names(cases) <- "season_cases_u5"
writeRaster(cases, "data/derived/season_cases_u5.tif", overwrite = TRUE)
writeRaster(season_inc, "data/derived/season_incidence_u5.tif", overwrite = TRUE)

# ---- district table -----------------------------------------------------------
ex <- extract(c(u5, cases), vect(adm2), fun = sum, na.rm = TRUE)
ds <- data.frame(region = adm2$region, district = adm2$shapeName, key = adm2$key,
                 children_3_59m = round(ex[, 2]), season_cases = round(ex[, 3])) |>
  left_join(dist |> transmute(key = adm2_key, season_inc_u5 = round(season_inc, 3),
                              annual_inc_u5 = round(annual_inc, 3), season_share_u5 = round(season_share_u5, 3)),
            by = "key") |>
  select(-key) |> arrange(desc(season_cases))
write.csv(ds, "data/derived/district_burden.csv", row.names = FALSE)

reg <- ds |> group_by(region) |>
  summarise(children = sum(children_3_59m), season_cases = sum(season_cases),
            per_child = round(season_cases / children, 3), .groups = "drop")
write.csv(reg, "data/derived/region_burden.csv", row.names = FALSE)
cat("By region (observed under-5 confirmed cases, Jun-Oct, DHIMS 2020-24):\n")
print(as.data.frame(reg), row.names = FALSE)
cat("Total children 3-59 m:", sum(reg$children), "| season cases:", sum(reg$season_cases), "\n")
cat("\nTop districts by season cases:\n"); print(head(ds, 8), row.names = FALSE)
cat("Highest per child:\n"); print(head(arrange(ds, desc(season_inc_u5)), 5), row.names = FALSE)
cat("Lowest per child:\n");  print(head(arrange(ds, season_inc_u5), 5), row.names = FALSE)

# ---- map ----------------------------------------------------------------------
png("out/burden_map.png", width = 2000, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 6))
plot(st_geometry(adm2), col = hcl.colors(9, "Reds", rev = TRUE)[cut(adm2$season_inc, 9)], border = "white", lwd = 0.4,
     main = "Under-5 confirmed cases per child, Jun-Oct (DHIMS 2020-24)")
plot(st_geometry(nr), add = TRUE, lwd = 1)
legend("bottomleft", levels(cut(adm2$season_inc, 9)), fill = hcl.colors(9, "Reds", rev = TRUE), bty = "n", cex = 0.6)
plot(log10(cases + 0.01), col = hcl.colors(50, "Reds", rev = TRUE), axes = FALSE, box = FALSE,
     main = "SMC-season under-5 cases per pixel (log10)")
plot(st_geometry(adm2), add = TRUE, border = "grey50", lwd = 0.3); plot(st_geometry(nr), add = TRUE)
dev.off()
cat("done\n")
