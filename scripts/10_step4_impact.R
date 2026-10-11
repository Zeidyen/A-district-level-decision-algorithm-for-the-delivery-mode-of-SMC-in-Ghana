# Step 4 + 5: cases averted by SMC under routine vs campaign delivery, and the
# green / amber / red map of where routine delivery could replace the campaign.
#
# Usage: Rscript scripts/10_step4_impact.R [scenario]   (default: main; see SCENARIOS)
#
# Per pixel, per SMC cycle (children 3-59 m, four cycles Jun-Sep):
#   campaign coverage   C_c, by region: Upper East 0.87 (0.867-0.895, Ayamba 2025);
#                       other regions 0.87 with a wide 0.75-0.95 (PLACEHOLDER until
#                       NMEP district data), since R scales as 1 / C_c
#   routine coverage    C_r = reach x P1_region x (1 - k x q)
#       P1_region : share ever reached by routine services = DHS 2022 Penta1
#                   (12-23 m, survey-weighted; checked against the data in the audit)
#       q         : per-contact probability that a monthly routine contact slips by
#                   >= 4 weeks or never happens (step 2 map, modelled directly)
#       k         : how much more often a routine malaria-product contact slips than
#                   a Penta contact. Lognormal, median 2, 95% ~1.4-2.8 (Ghana RTS,S
#                   vs Penta / MR dropout: Asante Lancet 2024; WHO SAGE report 2021)
#       reach     : 1 in the main analysis (scenario reach90: 0.9)
#   cases averted = baseline (no-SMC) season cases x PE x C
#   ratio R       = routine averted / campaign averted
# Expected cases averted are linear in per-cycle coverage (direct effect only), so
# whether the same children are missed each cycle does not change R.
#
# Baseline without SMC: step 3 / (1 - PE x C_c). It cancels from R; totals carry
# step 3's AGE_RATIO placeholder and are not interpretable until that is sourced.
#
# Classes (thresholds adopted 1 Oct 2026): green R >= 0.9, amber 0.7-0.9, red < 0.7,
# by the posterior mean of the district ratio, reported with P(class). 1,000 joint
# draws (step 2 surface, k, PE, C_c), in chunks to bound memory; seeded.
#
# Outputs (main): data/derived/district_step4.csv, step4_{ratio_mean,p_green,p_red,
# class}.tif, step4_k_stress.csv, out/step4_map.png. Other scenarios write only
# data/derived/district_step4_<scenario>.csv.

source("scripts/lib_inla.R")

SCENARIOS <- list(
  main    = list(outcome = "slip56", tt = "motor", region = FALSE, reach = 1, age = FALSE, adh_f = NA,
                 spatial = TRUE, covars = TRUE, ml = FALSE, prior_range = 30),
  thr49   = list(outcome = "slip49"),        # partial delays count: gap > 7 weeks
  thr42   = list(outcome = "slip42"),        # gap > 6 weeks
  walking = list(tt = "walk"),               # walking instead of motorised time
  region  = list(region = TRUE),             # region fixed effects in step 2
  reach90 = list(reach = 0.9),               # routine first contact reaches 90% of Penta1
  age     = list(age = TRUE),                # routine contact falls off at 24-59 months
  adh     = list(adh_f = 0.5),               # lower day-2/3 adherence under routine delivery
  adh0    = list(adh_f = 0),                 # ... and a partial course protects not at all
  age_adh = list(age = TRUE, adh_f = 0.5),   # both
  nospatial  = list(spatial = FALSE),                          # covariates, no spatial field
  regiononly = list(spatial = FALSE, covars = FALSE, region = TRUE), # regional means only
  stacked    = list(ml = TRUE, covars = FALSE),  # stacked ML covariate + spatial field (25/26)
  prior15    = list(prior_range = 15)            # field prior P(range < 15 km) = 0.05 (main: 30 km);
)                                                #   a prior sensitivity, not one of the 14 variants
# stacked: the step 2 covariates are replaced by the stacked machine-learning prediction
#   (25_ml_stack.py: out-of-fold for the clusters, learners on all clusters for the grid).
# age: routine contact at 24-59 m relative to 12-23 m, from vitamin A in the last 6 months
#   (DHS 2022, north, survey-weighted; 18_age_reach.R): 0.78 (0.73-0.84), drawn per draw.
#   SMC-eligible children split by WorldPop: 3-11 m = 0.75 x age 0; 12-23 m = 1/4 of ages
#   1-4; 24-59 m = 3/4 of ages 1-4. The door-to-door campaign is assumed to reach all ages.
# adh_f: day-2/3 adherence. Campaign 0.845 (Ayamba 2025, mean of 4 cycles). Routine scaled
#   by the Nigerian non-door-to-door vs door-to-door ratio 0.8674 / 0.9850 (Huang 2024,
#   Malar J 23:131) -> 0.744. A partial course is assumed to give a fraction adh_f of full
#   protection (ASSUMPTION: no sourced estimate; 0.5 main, 0 worst case). Routine coverage
#   is multiplied by (a_r + (1 - a_r) f) / (a_c + (1 - a_c) f).
args <- commandArgs(trailingOnly = TRUE)
SCEN <- if (length(args)) args[1] else "main"
stopifnot(SCEN %in% names(SCENARIOS))
sc <- modifyList(SCENARIOS$main, SCENARIOS[[SCEN]])
suffix <- if (SCEN == "main") "" else paste0("_", SCEN)
cat("scenario:", SCEN, "|", paste(names(sc), unlist(sc), sep = "=", collapse = " "), "\n")

N_DRAWS <- 1000; CHUNK <- 200
PE <- c(mean = 0.88, lo = 0.79, hi = 0.94)                     # Cairns 2021, days 0-28
CC_UE <- c(mean = 0.87, lo = 0.867, hi = 0.895)                # Ayamba 2025
CC_OTHER <- c(mean = 0.87, lo = 0.75, hi = 0.95)               # placeholder
K_MEDIAN <- 2; K_SDLOG <- (log(3) - log(1.5)) / 3.92
GREEN <- 0.9; RED <- 0.7
K_STRESS <- c(1, 1.5, 2, 3)
ADH_C <- 0.845; ADH_R <- ADH_C * 0.8674 / 0.9850
AGE_RATIO_OLD <- read.csv("data/derived/age_reach.csv")
P1 <- c(Northern = 0.895, Savannah = 0.962, `North East` = 0.955,
        `Upper East` = 0.982, `Upper West` = 0.973)

# Beta matching a mean and 95% interval (sd ~ width / 3.92).
rbeta_ci <- function(n, m, lo, hi) {
  v <- ((hi - lo) / 3.92)^2; a <- m * (1 - m) / v - 1
  rbeta(n, m * a, (1 - m) * a)
}

# ---- inputs -------------------------------------------------------------------
cl <- load_clusters(sc$outcome)
if (sc$ml) {
  mlf <- read.csv("data/derived/ml_oof_full.csv")
  cl$ml <- mlf$ml_stack[match(cl$cluster, mlf$cluster)]
  cl <- cl[!is.na(cl$ml), ]
}
m <- fit_model(cl, make_mesh(cl), list(tt = sc$tt, spatial = sc$spatial, region = sc$region, covars = sc$covars,
                                       ml = sc$ml, prior_range = sc$prior_range))
if (sc$spatial) cat(sprintf("field range %.0f km (95%% CrI %.0f-%.0f) | WAIC %.1f\n",
    m$fit$summary.hyperpar["Range for s", "mean"], m$fit$summary.hyperpar["Range for s", "0.025quant"],
    m$fit$summary.hyperpar["Range for s", "0.975quant"], m$fit$waic$waic))
cat("step 2 refitted; OR per doubling",
    if (sc$covars) round(exp(m$fit$summary.fixed["x1", "mean"]), 2) else "n/a (no covariates)",
    "| observed north q", round(sum(cl$y[cl$region %in% NORTH]) / sum(cl$n[cl$region %in% NORTH]), 3), "\n")
gd <- north_grid()
if (sc$ml) {
  mg <- read.csv("data/derived/ml_grid_pred.csv")
  gd$ml <- mg$ml[match(gd$cell, mg$cell)]
  gd$ml[is.na(gd$ml)] <- mean(gd$ml, na.rm = TRUE)            # 2 cells lack a covariate
}
ri <- match(gd$region, NORTH)

tt <- rast("data/derived/tt_motor_any_gha.tif")
rate <- rast("data/derived/season_incidence_u5.tif")          # step 3, per child per season
kids <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") * 0.75 +
        rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
kids <- resample(kids, tt, method = "bilinear"); kids[is.na(kids) | kids < 0] <- 0
gd$kids <- extract(kids, gd[, c("x", "y")])[, 2]
gd$rate <- extract(rate, gd[, c("x", "y")])[, 2]
gd$rate[is.na(gd$rate)] <- 0
obs <- gd$rate * gd$kids                                       # season cases under current SMC
k0 <- resample(rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif"), tt, method = "bilinear")
k1 <- resample(rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif"), tt, method = "bilinear")
a0 <- extract(k0, gd[, c("x", "y")])[, 2]; a1 <- extract(k1, gd[, c("x", "y")])[, 2]
a0[is.na(a0) | a0 < 0] <- 0; a1[is.na(a1) | a1 < 0] <- 0
share_old <- ifelse(0.75 * a0 + a1 > 0, 0.75 * a1 / (0.75 * a0 + a1), 0)   # 24-59 m among 3-59 m
adh_mult <- if (is.na(sc$adh_f)) 1 else
  (ADH_R + (1 - ADH_R) * sc$adh_f) / (ADH_C + (1 - ADH_C) * sc$adh_f)
cat(sprintf("age-weighting: share of children aged 24-59 m %.2f | adherence multiplier %.3f\n",
            sum(share_old * gd$kids) / sum(gd$kids), adh_mult))
w <- gd$kids

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)])
adm2 <- adm2[adm2$region %in% NORTH, ]
did <- sapply(st_intersects(st_as_sf(gd, coords = c("x", "y"), crs = 4326), adm2),
              function(i) if (length(i)) i[1] else NA)
inD <- !is.na(did); nD <- nrow(adm2)

# ---- draws ----------------------------------------------------------------------
set.seed(20261001)
kd <- exp(rnorm(N_DRAWS, log(K_MEDIAN), K_SDLOG))
pe <- rbeta_ci(N_DRAWS, PE["mean"], PE["lo"], PE["hi"])
ccm <- t(sapply(NORTH, function(r) {                          # region x draw
  p <- if (r == "Upper East") CC_UE else CC_OTHER
  rbeta_ci(N_DRAWS, p["mean"], p["lo"], p["hi"])
}))

ra <- if (sc$age) pmin(1, rnorm(N_DRAWS, AGE_RATIO_OLD$ratio, (AGE_RATIO_OLD$hi - AGE_RATIO_OLD$lo) / 3.92)) else rep(1, N_DRAWS)
sumR <- nG <- nRd <- sumCr <- sumq <- numeric(nrow(gd))
AC <- AR <- matrix(0, nD, N_DRAWS)
ARk <- setNames(lapply(K_STRESS, function(k) matrix(0, nD, N_DRAWS)), K_STRESS)
dsum <- function(M) {                                          # district sums per draw
  out <- matrix(0, nD, ncol(M)); r <- rowsum(M[inD, , drop = FALSE], did[inD])
  out[as.integer(rownames(r)), ] <- r; out
}

for (ch in seq_len(N_DRAWS / CHUNK)) {
  j <- (ch - 1) * CHUNK + seq_len(CHUNK)
  set.seed(20261001 + ch)                                      # iid draws inside draw_p
  q <- draw_p(m, gd, CHUNK, seed = 20261001L + ch)             # pixels x draws
  ccp <- ccm[ri, j]                                            # campaign coverage per pixel
  qk <- sweep(q, 2, kd[j], "*"); qk[qk > 1] <- 1               # (pmin would drop dims)
  agef <- 1 - outer(share_old, 1 - ra[j])                      # pixels x draws
  Cr <- sc$reach * P1[gd$region] * (1 - qk) * agef * adh_mult
  R <- Cr / ccp
  sumR <- sumR + rowSums(R); nG <- nG + rowSums(R >= GREEN); nRd <- nRd + rowSums(R < RED)
  sumCr <- sumCr + rowSums(Cr); sumq <- sumq + rowSums(q)
  base <- obs / (1 - sweep(ccp, 2, pe[j], "*"))                # no-SMC baseline
  AC[, j] <- dsum(base * sweep(ccp, 2, pe[j], "*"))
  AR[, j] <- dsum(base * sweep(Cr, 2, pe[j], "*"))
  for (k in K_STRESS) {                                        # k fixed, same draws
    qf <- k * q; qf[qf > 1] <- 1
    ARk[[as.character(k)]][, j] <- dsum(base * sweep(sc$reach * P1[gd$region] * (1 - qf), 2, pe[j], "*"))
  }
  rm(q, ccp, qk, Cr, R, base); invisible(gc())
  cat("  draws", max(j), "of", N_DRAWS, "\n")
}

# ---- districts ----------------------------------------------------------------
cls <- function(r) c("red", "amber", "green")[ifelse(r >= GREEN, 3L, ifelse(r >= RED, 2L, 1L))]
RD <- AR / AC
kidsD <- tapply(w[inD], factor(did[inD], seq_len(nD)), sum)
crD <- tapply((sumCr / N_DRAWS * w)[inD], factor(did[inD], seq_len(nD)), sum) / kidsD
ds <- data.frame(region = adm2$region, district = adm2$shapeName,
                 children_3_59m = round(kidsD), routine_cov = round(crD, 3),
                 ratio = round(rowMeans(RD), 3),
                 lo = round(apply(RD, 1, quantile, 0.025), 3), hi = round(apply(RD, 1, quantile, 0.975), 3),
                 p_green = round(rowMeans(RD >= GREEN), 2), p_amber = round(rowMeans(RD >= RED & RD < GREEN), 2),
                 p_red = round(rowMeans(RD < RED), 2),
                 class = cls(rowMeans(RD)),
                 averted_campaign = round(rowMeans(AC)), averted_routine = round(rowMeans(AR)),
                 shortfall = round(rowMeans(AC - AR)))
ds <- ds[order(ds$ratio), ]
write.csv(ds, sprintf("data/derived/district_step4%s.csv", suffix), row.names = FALSE)

tab <- table(factor(ds$class, c("green", "amber", "red")))
conf <- c(green = sum(ds$p_green >= 0.8), red = sum(ds$p_red >= 0.8))
cat("\nDistricts by class:", paste(names(tab), tab, collapse = " | "),
    "\n  with P(class) >= 0.8: green", conf["green"], "| red", conf["red"], "\n")
cat("Children 3-59 m by district class:",
    paste(names(tab), round(tapply(ds$children_3_59m, factor(ds$class, names(tab)), sum)), collapse = " | "), "\n")
cat("Routine per-cycle coverage, child-weighted:",
    paste(NORTH, round(sapply(NORTH, function(r) sum((sumCr / N_DRAWS * w)[gd$region == r]) / sum(w[gd$region == r])), 2),
          collapse = " | "), "| all", round(sum(sumCr / N_DRAWS * w) / sum(w), 2), "\n")
cat(sprintf("Routine shortfall: %.0f%% of campaign cases averted; %.0f%% of it in red districts\n",
            100 * sum(ds$shortfall) / sum(ds$averted_campaign),
            100 * sum(ds$shortfall[ds$class == "red"]) / sum(ds$shortfall)))
cat("\nLowest ratios:\n"); print(head(ds[, c("region", "district", "routine_cov", "ratio", "lo", "hi", "p_green", "p_red", "class")], 8), row.names = FALSE)
cat("Highest:\n");   print(tail(ds[, c("region", "district", "routine_cov", "ratio", "lo", "hi", "p_green", "p_red", "class")], 5), row.names = FALSE)

if (SCEN != "main") quit(save = "no")

# ---- stress test (main only): k fixed --------------------------------------------
st <- do.call(rbind, lapply(K_STRESS, function(k) {
  rk <- rowMeans(ARk[[as.character(k)]] / AC)
  ck <- cls(rk)
  data.frame(k = k, green = sum(ck == "green"), amber = sum(ck == "amber"), red = sum(ck == "red"),
             pct_children_red = round(100 * sum(kidsD[ck == "red"]) / sum(kidsD)))
}))
write.csv(st, "data/derived/step4_k_stress.csv", row.names = FALSE)
cat("\nStress test (k fixed):\n"); print(st, row.names = FALSE)

# Region break-even k: k at which the region's child- and burden-weighted ratio
# (posterior-mean q, campaign coverage at its mean) falls to the green line.
qbar <- sumq / N_DRAWS
be <- sapply(NORTH, function(r) {
  i <- gd$region == r
  f <- function(k) sum(sc$reach * P1[r] * (1 - pmin(1, k * qbar[i])) * obs[i]) / sum(obs[i]) /
                   (if (r == "Upper East") CC_UE else CC_OTHER)["mean"] - GREEN
  if (f(0) < 0) 0 else if (f(20) > 0) Inf else uniroot(f, c(0, 20))$root
})
cat("Break-even k by region:", paste(names(be), round(be, 2), collapse = " | "), "\n")
write.csv(data.frame(region = names(be), breakeven_k = round(be, 2)), "data/derived/step4_breakeven.csv",
          row.names = FALSE)

# ---- pixel outputs + map --------------------------------------------------------
out <- rast(tt); out[] <- NA
pix <- cellFromXY(out, as.matrix(gd[, c("x", "y")]))        # gd$cell indexes the cropped grid
put <- function(v, nm) { r <- out; r[pix] <- v; names(r) <- nm; r }
pc <- ifelse(sumR / N_DRAWS >= GREEN, 3L, ifelse(sumR / N_DRAWS >= RED, 2L, 1L))
rr <- c(put(sumR / N_DRAWS, "ratio_mean"), put(nG / N_DRAWS, "p_green"), put(nRd / N_DRAWS, "p_red"),
        put(pc, "class"))
nr <- adm1[sub(" Region$", "", adm1$shapeName) %in% NORTH, ]
rr <- crop(rr, ext(vect(nr)))
for (i in 1:4) writeRaster(rr[[i]], sprintf("data/derived/step4_%s.tif", names(rr)[i]), overwrite = TRUE)

cols3 <- c("#c0392b", "#f0ad4e", "#2e8b57")
png("out/step4_map.png", width = 2400, height = 1000, res = 150)
par(mfrow = c(1, 3), mar = c(1, 1, 3, 6))
plot(rr[["class"]], col = cols3, type = "classes", levels = c("red", "amber", "green"),
     main = "Pixel class (mean ratio)", axes = FALSE, box = FALSE)
plot(st_geometry(adm2), add = TRUE, border = "grey40", lwd = 0.3); plot(st_geometry(nr), add = TRUE)
plot(rr[["p_green"]], col = hcl.colors(50, "Greens", rev = TRUE), range = c(0, 1),
     main = "P(routine >= 90% of campaign)", axes = FALSE, box = FALSE)
plot(st_geometry(adm2), add = TRUE, border = "grey40", lwd = 0.3); plot(st_geometry(nr), add = TRUE)
adm2$class <- ds$class[match(adm2$shapeName, ds$district)]
plot(st_geometry(adm2), col = cols3[match(adm2$class, c("red", "amber", "green"))],
     border = "white", main = "District class")
plot(st_geometry(nr), add = TRUE, lwd = 1.5)
legend("bottomleft", c("green: >= 90%", "amber: 70-90%", "red: < 70%"), fill = rev(cols3), bty = "n", cex = 0.9)
dev.off()
cat("done\n")
