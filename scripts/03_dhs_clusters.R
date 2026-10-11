# Step 2 input: cluster-level routine-contact indicators from Ghana DHS 2022,
# joined to cluster GPS and to the step 1 travel-time surfaces.
#
# Indicators (children alive, standard DHS definitions):
#   penta1, penta3  : received DPT/Penta dose 1 / 3, children 12-23 months
#                     (card date, card mark, or mother's report all count)
#   vita6m          : vitamin A in the last 6 months, children 6-59 months
# Dropout Penta1->3 is derived from penta1 and penta3.
#
# GPS is displaced (urban <=2 km, rural <=5 km, 1% <=10 km), restricted to the
# district. Travel time is therefore averaged over a buffer around each cluster
# (2 km urban, 5 km rural), not read at the point, as DHS guidance recommends.
#
# Output: data/derived/dhs_clusters.csv (all Ghana clusters; filter for the
# north at model time, since southern clusters help estimate covariate effects).

suppressMessages({library(haven); library(sf); library(terra); library(dplyr)})
sf_use_s2(FALSE)

kr <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA")
got <- function(x) as.integer(x %in% 1:3)    # 1 date on card, 2 reported, 3 marked
kr <- kr |>
  filter(b5 == 1) |>
  mutate(age = hw1, age_m = v008 - b3,       # hw1 is missing for some; use CMC age
         age_m = ifelse(is.na(age_m), age, age_m),
         penta1 = got(h3), penta3 = got(h7), vita6m = as.integer(h34 == 1),
         w = v005 / 1e6)

clus <- kr |>
  group_by(cluster = v001) |>
  summarise(
    n_12_23   = sum(age_m >= 12 & age_m <= 23),
    penta1    = sum(penta1[age_m >= 12 & age_m <= 23]),
    penta3    = sum(penta3[age_m >= 12 & age_m <= 23]),
    n_6_59    = sum(age_m >= 6 & age_m <= 59 & !is.na(h34)),
    vita6m    = sum(vita6m[age_m >= 6 & age_m <= 59 & !is.na(h34)]),
    wealth    = mean(v190),                    # 1 poorest .. 5 richest
    mother_sec = mean(v106 >= 2),              # secondary or higher
    weight    = mean(w),
    .groups = "drop")

ge <- st_read("data/raw/dhs/GHGE8AFL/GHGE8AFL.shp", quiet = TRUE) |>
  select(cluster = DHSCLUST, region = ADM1NAME, urban = URBAN_RURA, LATNUM, LONGNUM)
d <- inner_join(ge, clus, by = "cluster")

# Covariates averaged over the displacement buffer, weighted by under-5
# population, from the Ghana-wide rasters of 05_access_ghana.R (run that first).
# The weighting matters: a plain buffer mean counts empty land as much as the
# villages where the children live, and overstated northern cluster travel times
# by about half (Northern median 51 vs 33 min), which made the step 2 map
# under-predict. Falls back to the plain mean if a buffer holds no population.
buf <- vect(st_buffer(d, ifelse(d$urban == "U", 2000, 5000) / 111320))
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") +
      rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
u5 <- resample(u5, rast("data/derived/tt_walk_any_gha.tif"))
u5[is.na(u5) | u5 < 0] <- 0
for (nm in c("tt_walk_any", "tt_motor_any", "logpopdens")) {
  r <- rast(sprintf("data/derived/%s_gha.tif", nm))
  e <- extract(c(r * u5, u5), buf, fun = sum, na.rm = TRUE)
  plain <- extract(r, buf, fun = mean, na.rm = TRUE)[, 2]
  d[[nm]] <- ifelse(e[, 3] > 0, e[, 2] / e[, 3], plain)
}
cat("clusters with no population in buffer (plain mean used):", sum(e[, 3] == 0), "\n")

out <- st_drop_geometry(d)
write.csv(out, "data/derived/dhs_clusters.csv", row.names = FALSE)

# Weighted regional check (children-level, survey weights).
north <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
cat("\nWeighted coverage by region (children 12-23 m for Penta; 6-59 m for vit A):\n")
reg <- kr |> mutate(region = as.character(as_factor(v024))) |>
  group_by(region) |>
  summarise(
    penta1 = weighted.mean(penta1[age_m %in% 12:23], w[age_m %in% 12:23]),
    penta3 = weighted.mean(penta3[age_m %in% 12:23], w[age_m %in% 12:23]),
    dropout = 1 - penta3 / penta1,
    vita6m = weighted.mean(vita6m[age_m %in% 6:59 & !is.na(h34)], w[age_m %in% 6:59 & !is.na(h34)]),
    n_12_23 = sum(age_m %in% 12:23), .groups = "drop") |>
  mutate(across(penta1:vita6m, ~ round(100 * .x, 1)))
print(as.data.frame(reg), row.names = FALSE)

n <- out[out$region %in% north, ]
cat("\nNorthern clusters:", nrow(n), "| children 12-23 m:", sum(n$n_12_23),
    "| clusters with 0 children 12-23 m:", sum(n$n_12_23 == 0), "\n")
cat("Spearman, cluster Penta3 share vs walking time (north):",
    round(cor(n$penta3 / n$n_12_23, n$tt_walk_any, method = "spearman", use = "complete"), 2), "\n")
cat("Spearman, cluster vit A share vs walking time (north):",
    round(cor(n$vita6m / n$n_6_59, n$tt_walk_any, method = "spearman", use = "complete"), 2), "\n")

# Other routine-contact indicators vs walking time (northern clusters), so every
# correlation quoted in the design note and summary comes from code:
# MCV2 (24-35 m), fever with medical treatment, mother reports distance as a big
# problem (v467d == 1, one per woman), and household wealth quintile.
kr2 <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA") |> filter(b5 == 1) |> mutate(age_m = v008 - b3)
ex <- kr2 |> group_by(cluster = v001) |> summarise(
  mcv2 = mean(h9a[age_m %in% 24:35] %in% 1:3), n_mcv2 = sum(age_m %in% 24:35),
  fever_tx = mean(h32z[h22 %in% 1 & !is.na(h32z)] == 1),
  dist_problem = mean(v467d[!duplicated(caseid)] == 1, na.rm = TRUE), .groups = "drop")
nn <- inner_join(out[out$region %in% north, ], ex, by = "cluster")
sp <- function(v) round(cor(v, nn$tt_walk_any, method = "spearman", use = "complete"), 2)
cat("Spearman with walking time (north):",
    "MCV2", sp(ifelse(nn$n_mcv2 > 0, nn$mcv2, NA)), "| fever treated", sp(nn$fever_tx),
    "| distance a big problem", sp(nn$dist_problem), "| wealth", sp(nn$wealth), "\n")
cat("Motorised time at the 5-min floor:", round(100 * mean(out$tt_motor_any < 5)), "% of clusters\n")
