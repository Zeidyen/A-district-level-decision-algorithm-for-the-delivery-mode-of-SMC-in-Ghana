# Does rainfall make routine contacts slip, and more so for remote children?
#
# Unit: one inter-dose gap (Penta 1->2 and 2->3), northern children 12-35 m with
# usable card dates (same rules as 04_timeliness.R).
#   slip  = the gap was longer than 8 weeks (the next monthly contact was missed
#           by a month or more)
#   rain  = CHIRPS rainfall (mm) at the cluster, averaged over the displacement
#           buffer, in the calendar month the next dose fell due (previous dose
#           + 28 days)
#
# Model: slip ~ rain x log2(walking time) + gap number + year + region
#              + (1 | cluster) + (1 | child)
# Year enters as a fixed effect because 2020 (COVID-19) disrupted services.
# Rainfall and the farming season coincide, so "rain" here means the wet season
# as a whole, not rain alone.

suppressMessages({library(haven); library(dplyr); library(tidyr); library(sf);
                  library(terra); library(glmmTMB)})
sf_use_s2(FALSE)

north <- c("northern", "savannah", "north east", "upper east", "upper west")
day <- function(d, m, y) {
  ok <- !is.na(d) & d >= 1 & d <= 31 & m >= 1 & m <= 12 & y >= 2015 & y <= 2023
  out <- rep(as.Date(NA), length(d))
  out[ok] <- as.Date(sprintf("%04d-%02d-%02d", as.integer(y[ok]), as.integer(m[ok]),
                             as.integer(d[ok])), optional = TRUE)
  out
}
kr <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA") |>
  mutate(age_m = v008 - b3, region = as.character(as_factor(v024))) |>
  filter(b5 == 1, age_m >= 12, age_m <= 35, region %in% north) |>
  mutate(child = row_number(), dob = day(b17, b1, b2),
         p1 = as.Date(ifelse(h3 == 1, day(h3d, h3m, h3y), NA)),
         p2 = as.Date(ifelse(h5 == 1, day(h5d, h5m, h5y), NA)),
         p3 = as.Date(ifelse(h7 == 1, day(h7d, h7m, h7y), NA)),
         ok = !is.na(dob) & !is.na(p1) & !is.na(p2) & !is.na(p3) &
              as.numeric(p3 - dob) >= 70 & p2 > p1 & p3 > p2) |>
  filter(ok)

gaps <- bind_rows(
  transmute(kr, child, cluster = v001, region, gap = "1->2", prev = p1, nxt = p2),
  transmute(kr, child, cluster = v001, region, gap = "2->3", prev = p2, nxt = p3)) |>
  mutate(slip = as.integer(as.numeric(nxt - prev) > 56),
         due = prev + 28, due_m = format(due, "%Y-%m"), year = format(due, "%Y"),
         due_month = as.integer(format(due, "%m")))

# Rainfall in the due month at each cluster (buffer mean, as for the covariates).
cl <- read.csv("data/derived/dhs_clusters.csv")
cl <- cl[cl$cluster %in% gaps$cluster, ]
buf <- st_buffer(st_as_sf(cl, coords = c("LONGNUM", "LATNUM"), crs = 4326),
                 ifelse(cl$urban == "U", 2000, 5000) / 111320)
rain <- rast("data/derived/chirps_monthly_gha.tif")
rm_ <- extract(rain, vect(buf), fun = mean, na.rm = TRUE)[, -1]
rm_ <- cbind(cluster = cl$cluster, rm_) |>
  pivot_longer(-cluster, names_to = "due_m", values_to = "rain")
gaps <- gaps |> left_join(rm_, by = c("cluster", "due_m")) |>
  left_join(cl[, c("cluster", "tt_walk_any")], by = "cluster") |>
  mutate(log2_tt = log2(pmax(tt_walk_any, 5)),
         log2_tt_c = log2_tt - log2(60),     # centred at 60 min walking
         rain100 = rain / 100,
         wet = due_month %in% 6:9)
cat(nrow(gaps), "gaps from", n_distinct(gaps$child), "children;",
    sum(is.na(gaps$rain)), "without rainfall\n\n")

# Descriptive: slip by season and distance.
gaps$far <- ifelse(gaps$tt_walk_any >= 60, ">=60 min walk", "<60 min walk")
cat("Share of gaps that slipped (%), dose due in Jun-Sep vs other months:\n")
print(gaps |> group_by(far, season = ifelse(wet, "Jun-Sep", "other")) |>
        summarise(n = n(), slip = round(100 * mean(slip)), .groups = "drop") |>
        pivot_wider(names_from = season, values_from = c(n, slip)) |> as.data.frame(),
      row.names = FALSE)

cat("\nBy rainfall in the due month:\n")
gaps$rain_band <- cut(gaps$rain, c(-Inf, 25, 100, 200, Inf), c("<25 mm", "25-100", "100-200", ">200"))
print(gaps |> group_by(rain_band) |> summarise(n = n(), slip = round(100 * mean(slip))) |>
        as.data.frame(), row.names = FALSE)

# Models.
ors <- function(m, terms) {
  b <- fixef(m)$cond; V <- vcov(m)$cond
  for (t in terms) cat(sprintf("  %-22s %.2f (%.2f-%.2f)\n", t, exp(b[t]),
                               exp(b[t] - 1.96 * sqrt(V[t, t])), exp(b[t] + 1.96 * sqrt(V[t, t]))))
}
m1 <- glmmTMB(slip ~ rain100 + log2_tt_c + gap + year + region + (1 | cluster) + (1 | child),
              data = gaps, family = binomial)
cat("\nMain effects (OR):\n"); ors(m1, c("rain100", "log2_tt_c"))
b1 <- fixef(m1)$cond; V1 <- vcov(m1)$cond
write.csv(data.frame(term = c("rain100", "log2_tt_c"), or = exp(b1[c("rain100", "log2_tt_c")]),
                     lo = exp(b1[c("rain100", "log2_tt_c")] - 1.96 * sqrt(diag(V1)[c("rain100", "log2_tt_c")])),
                     hi = exp(b1[c("rain100", "log2_tt_c")] + 1.96 * sqrt(diag(V1)[c("rain100", "log2_tt_c")]))),
          "data/derived/rain_or.csv", row.names = FALSE)
m2 <- glmmTMB(slip ~ rain100 * log2_tt_c + gap + year + region + (1 | cluster) + (1 | child),
              data = gaps, family = binomial)
cat("\nWith interaction (rain OR at 60 min walk; interaction = change in rain OR per doubling):\n")
ors(m2, c("rain100", "log2_tt_c", "rain100:log2_tt_c"))
cat(sprintf("  AIC: main %.1f | interaction %.1f\n", AIC(m1), AIC(m2)))

# Rain OR (per 100 mm) at chosen walking times, from the interaction model.
b <- fixef(m2)$cond; V <- vcov(m2)$cond
cat("\nRain OR per 100 mm, at walking time of:\n")
for (t in c(15, 30, 60, 120, 240)) {
  x <- log2(t) - log2(60); g <- c(1, x); nm <- c("rain100", "rain100:log2_tt_c")
  est <- sum(g * b[nm]); se <- sqrt(as.numeric(t(g) %*% V[nm, nm] %*% g))
  cat(sprintf("  %3d min: %.2f (%.2f-%.2f)\n", t, exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se)))
}

# Season as a binary, which is what SMC timing is defined by.
m3 <- glmmTMB(slip ~ wet * log2_tt_c + gap + year + region + (1 | cluster) + (1 | child),
              data = gaps, family = binomial)
cat("\nJun-Sep vs other months:\n"); ors(m3, c("wetTRUE", "wetTRUE:log2_tt_c"))

# Survivorship check. The gaps above come only from children who received all
# three doses. If rain makes children drop out rather than come late, that is
# invisible above. Children with a card-dated Penta1: did those whose Penta2 fell
# due in Jun-Sep more often never complete Penta3?
d1 <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA") |>
  mutate(age_m = v008 - b3, region = as.character(as_factor(v024))) |>
  filter(b5 == 1, age_m >= 12, age_m <= 35, region %in% north, h3 == 1) |>
  mutate(p1 = day(h3d, h3m, h3y), dropout = as.integer(h7 == 0),
         wet = as.integer(format(p1 + 28, "%m")) %in% 6:9) |>
  filter(!is.na(p1)) |>
  left_join(cl[, c("cluster", "tt_walk_any")], by = c("v001" = "cluster")) |>
  mutate(far = ifelse(tt_walk_any >= 60, ">=60 min walk", "<60 min walk"))
cat("\nNever received Penta3 (%), by when Penta2 fell due:\n")
print(d1 |> group_by(far, season = ifelse(wet, "Jun-Sep", "other")) |>
        summarise(n = n(), dropout = round(100 * mean(dropout), 1), .groups = "drop") |>
        as.data.frame(), row.names = FALSE)
m4 <- glmmTMB(dropout ~ wet * I(log2(pmax(tt_walk_any, 5)) - log2(60)) + region + (1 | v001),
              data = d1, family = binomial)
cat("Dropout model:\n"); ors(m4, names(fixef(m4)$cond)[2:3])
