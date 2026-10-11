# Break-even k by district, and a global sensitivity analysis of the routine : campaign ratio.
#
# Why: the class counts in step 4 rest on k (how much more often a routine SMC contact
# slips than a Penta contact), which no study measures for SMC. Break-even k turns the
# question round: how reliable would routine SMC have to be for routine delivery to reach
# 90% of the campaign's cases averted? A break-even below 1 means routine SMC would have to
# be MORE reliable than routine vaccination, which no evidence suggests. The global
# sensitivity analysis then asks which assumption drives the ratio in each region.
#
# Same model, data and draws as 10_step4_impact.R (main): the 1,000 step-2 draws use the
# same seeds. Within a district, campaign coverage and the no-SMC baseline are constant
# (region-level), so drug protection and the baseline cancel and
#   R_d = reach x P1 x m_adh / C_c x sum_i obs_i (1 - min(1, k q_i)) (1 - s_i (1 - ra)) / sum_i obs_i
# with obs = observed season cases (step 3 weights), s = share of children aged 24-59 m,
# ra = routine reach at 24-59 m relative to 12-23 m, m_adh = adherence multiplier.
# District sums are precomputed per draw on a grid of k, so everything below is exact up to
# linear interpolation in k (grid step 0.05).
#
# Break-even k, two settings:
#   best case  : no age gap (ra = 1), no adherence loss, reach 1 (as step4_breakeven.csv)
#   realistic  : age gap ra ~ N(0.78, CI 0.73-0.84) and adherence (partial course f = 0.5)
# Uncertainty from step-2 draws and campaign coverage; reported as median and 95% interval.
#
# Global sensitivity (N = 20,000; each step-2 draw used 20 times). Inputs, independent:
#   k ~ lognormal(median 2, 95% 1.4-2.8) | C_c as step 4 | ra ~ U(0.73, 1) (measured gap to
#   none) | routine adherence a_r ~ U(0.744, 0.845) (Nigerian gap to none) | partial-course
#   protection f ~ U(0, 1) | reach ~ U(0.9, 1) | step-2 surface (summarised by the region's
#   burden-weighted mean q in that draw). Drug protection (PE) cancels from the ratio.
# First-order indices: Var(E[Y | X]) / Var(Y), X binned into 20 quantile bins.
#
# Favourable-to-routine check: k = 1, no age gap, no adherence loss, reach 1, C_c at its mean
# and at its low end (0.75 outside Upper East).
#
# Outputs: data/derived/district_breakeven.csv, region_breakeven.csv, gsa_indices.csv,
#          gsa_region_ratio.csv, favourable_check.csv; out/breakeven_map.png (best case, 90% and
#          70% lines; the realistic setting is in the CSVs: every district is below 1),
#          out/gsa_chart.png

source("scripts/lib_inla.R")
N_DRAWS <- 1000; CHUNK <- 200
CC_UE <- c(mean = 0.87, lo = 0.867, hi = 0.895); CC_OTHER <- c(mean = 0.87, lo = 0.75, hi = 0.95)
# C_c: survey-measured where available (Upper East, Ayamba 2025), assumed elsewhere. Must agree with
# MEASURED_CC in 31_delivery_algorithm.R, whose 80% programme target is a separate quality gate.
K_MEDIAN <- 2; K_SDLOG <- (log(3) - log(1.5)) / 3.92
GREEN <- 0.9; RED <- 0.7
ADH_C <- 0.845; ADH_R <- ADH_C * 0.8674 / 0.9850
P1 <- c(Northern = 0.895, Savannah = 0.962, `North East` = 0.955, `Upper East` = 0.982, `Upper West` = 0.973)
AGE <- read.csv("data/derived/age_reach.csv")
KG <- seq(0, 6, by = 0.05); nK <- length(KG)
rbeta_ci <- function(n, m, lo, hi) {
  v <- ((hi - lo) / 3.92)^2; a <- m * (1 - m) / v - 1
  rbeta(n, m * a, (1 - m) * a)
}
ccpar <- function(r) if (r == "Upper East") CC_UE else CC_OTHER

# ---- inputs (as 10_step4_impact.R) -----------------------------------------------
cl <- load_clusters("slip56")
m <- fit_model(cl, make_mesh(cl), list(tt = "motor", spatial = TRUE, region = FALSE, covars = TRUE))
gd <- north_grid()
tt <- rast("data/derived/tt_motor_any_gha.tif")
rate <- rast("data/derived/season_incidence_u5.tif")
r0 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif"); r1 <- rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
ex <- function(r) { v <- extract(resample(r, tt, method = "bilinear"), gd[, c("x", "y")])[, 2]; v[is.na(v) | v < 0] <- 0; v }
a0 <- ex(r0); a1 <- ex(r1)
kids <- 0.75 * a0 + a1
rt <- extract(rate, gd[, c("x", "y")])[, 2]; rt[is.na(rt)] <- 0
obs <- rt * kids
s_old <- ifelse(kids > 0, 0.75 * a1 / kids, 0)

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)])
adm2 <- adm2[adm2$region %in% NORTH, ]
did <- sapply(st_intersects(st_as_sf(gd, coords = c("x", "y"), crs = 4326), adm2),
              function(i) if (length(i)) i[1] else NA)
inD <- !is.na(did) & obs > 0; nD <- nrow(adm2)
dreg <- adm2$region
fD <- factor(did[inD], seq_len(nD))
A  <- as.vector(tapply(obs[inD], fD, sum))                   # sum obs
Cs <- as.vector(tapply((obs * s_old)[inD], fD, sum))         # sum obs x s

# ---- district sums per draw on the k grid ------------------------------------------
S  <- array(0, c(nD, nK, N_DRAWS))                           # sum obs (1 - min(1, k q))
Sd <- array(0, c(nD, nK, N_DRAWS))                           # sum obs s (1 - min(1, k q))
QB <- matrix(0, nD, N_DRAWS)                                 # sum obs q
oi <- obs[inD]; osi <- (obs * s_old)[inD]; di <- did[inD]
for (ch in seq_len(N_DRAWS / CHUNK)) {
  j <- (ch - 1) * CHUNK + seq_len(CHUNK)
  set.seed(20261001 + ch)
  q <- draw_p(m, gd, CHUNK, seed = 20261001L + ch)[inD, , drop = FALSE]
  QB[, j] <- rowsum(q * oi, di, reorder = TRUE)[as.character(seq_len(nD)), ]
  for (ik in seq_len(nK)) {
    keep <- 1 - pmin(1, KG[ik] * q); dim(keep) <- dim(q)
    S[, ik, j]  <- rowsum(keep * oi, di, reorder = TRUE)[as.character(seq_len(nD)), ]
    Sd[, ik, j] <- rowsum(keep * osi, di, reorder = TRUE)[as.character(seq_len(nD)), ]
  }
  rm(q); invisible(gc()); cat("  draws", max(j), "of", N_DRAWS, "\n")
}

# Ratio for districts (or a region: pass summed rows) at given settings, one draw index j.
# Returns a vector over the k grid.
ratio_k <- function(Sx, Sdx, Ax, Cx, p1, cc, ra, madh, reach)
  reach * p1 * madh / cc * (Sx - (1 - ra) * Sdx) / Ax
interp_k <- function(rk, k) approx(KG, rk, xout = pmin(k, max(KG)), rule = 2)$y
# Smallest k at which the ratio falls to threshold T (ratio decreases in k).
breakeven <- function(rk, T) {
  if (rk[1] < T) return(0)
  if (rk[nK] >= T) return(Inf)
  i <- which(rk < T)[1]
  KG[i - 1] + (KG[i] - KG[i - 1]) * (rk[i - 1] - T) / (rk[i - 1] - rk[i])
}

set.seed(20261021)
madh_f <- function(ar, f) (ar + (1 - ar) * f) / (ADH_C + (1 - ADH_C) * f)
cc_d <- sapply(seq_len(nD), function(d) { p <- ccpar(dreg[d]); rbeta_ci(N_DRAWS, p["mean"], p["lo"], p["hi"]) })  # draw x district
ra_d <- pmin(1, rnorm(N_DRAWS, AGE$ratio, (AGE$hi - AGE$lo) / 3.92))
madh_real <- madh_f(ADH_R, 0.5)

# ---- break-even k by district ----------------------------------------------------
be <- function(setting, T) sapply(seq_len(nD), function(d) sapply(seq_len(N_DRAWS), function(j) {
  ra <- if (setting == "best") 1 else ra_d[j]; ma <- if (setting == "best") 1 else madh_real
  breakeven(ratio_k(S[d, , j], Sd[d, , j], A[d], Cs[d], P1[dreg[d]], cc_d[j, d], ra, ma, 1), T)
}))                                                          # draw x district
BEb <- be("best", GREEN); BEr <- be("real", GREEN); BEb7 <- be("best", RED)
q3 <- function(M, p) round(apply(M, 2, quantile, p), 2)
kev <- exp(rnorm(N_DRAWS, log(K_MEDIAN), K_SDLOG))
dist <- data.frame(region = dreg, district = adm2$shapeName,
                   children_3_59m = round(as.vector(tapply(kids[inD], fD, sum))),
                   mean_q = round(rowMeans(QB) / A, 3),
                   be_best = q3(BEb, 0.5), be_best_lo = q3(BEb, 0.025), be_best_hi = q3(BEb, 0.975),
                   be_real = q3(BEr, 0.5), be_real_lo = q3(BEr, 0.025), be_real_hi = q3(BEr, 0.975),
                   be_best_red = q3(BEb7, 0.5),
                   p_best_below1 = round(colMeans(BEb < 1), 2),
                   # For the delivery-mode algorithm (31): probability that break-even k lies below /
                   # above the evidence range for k (1.4-2.8, the lognormal's 95% range).
                   p_best_below_klo = round(colMeans(BEb < 1.4), 3), p_best_above_khi = round(colMeans(BEb > 2.8), 3),
                   p_real_below_klo = round(colMeans(BEr < 1.4), 3), p_real_above_khi = round(colMeans(BEr > 2.8), 3),
                   p_best_below_k = round(colMeans(BEb < kev), 2))   # P(green line missed | k as evidence)
dist <- dist[order(dist$be_best), ]
write.csv(dist, "data/derived/district_breakeven.csv", row.names = FALSE)

# Region break-even: posterior-mean sums, campaign coverage at its mean (as step4_breakeven.csv).
regbe <- sapply(NORTH, function(r) {
  i <- dreg == r
  rk <- ratio_k(rowMeans(colSums(S[i, , , drop = FALSE]), dims = 1), rowMeans(colSums(Sd[i, , , drop = FALSE]), dims = 1),
                sum(A[i]), sum(Cs[i]), P1[r], ccpar(r)["mean"], 1, 1, 1)
  rk2 <- ratio_k(rowMeans(colSums(S[i, , , drop = FALSE]), dims = 1), rowMeans(colSums(Sd[i, , , drop = FALSE]), dims = 1),
                 sum(A[i]), sum(Cs[i]), P1[r], ccpar(r)["mean"], AGE$ratio, madh_real, 1)
  c(best = breakeven(rk, GREEN), realistic = breakeven(rk2, GREEN), best_red = breakeven(rk, RED))
})
cat("\nBreak-even k by region (green line; best case | realistic; red line, best case):\n")
print(round(t(regbe), 2))
old <- read.csv("data/derived/step4_breakeven.csv")
cat("check vs step4_breakeven.csv (best case):", paste(old$region, old$breakeven_k, collapse = " | "), "\n")
write.csv(data.frame(region = NORTH, round(t(regbe), 2)), "data/derived/region_breakeven.csv", row.names = FALSE)

# Validation: reproduce the main step 4 district ratio (k, C_c, draws as main are not the same
# RNG stream, so compare posterior means at k = 2 and C_c mean rather than draw by draw).
st4 <- read.csv("data/derived/district_step4.csv")
chk <- sapply(seq_len(nD), function(d) mean(sapply(seq_len(N_DRAWS), function(j)
  interp_k(ratio_k(S[d, , j], Sd[d, , j], A[d], Cs[d], P1[dreg[d]], ccpar(dreg[d])["mean"], 1, 1, 1), 2))))
cat("validation: district ratio at k = 2 vs step 4 main posterior mean, correlation",
    round(cor(chk, st4$ratio[match(adm2$shapeName, st4$district)]), 3),
    "| mean abs diff", round(mean(abs(chk - st4$ratio[match(adm2$shapeName, st4$district)])), 3), "\n")

cat("\nDistricts by best-case break-even k (green line):\n")
cat("  below 1 (routine SMC would have to beat vaccination):", sum(dist$be_best < 1),
    "districts,", round(100 * sum(dist$children_3_59m[dist$be_best < 1]) / sum(dist$children_3_59m)), "% of children\n")
cat("  1 to 1.4:", sum(dist$be_best >= 1 & dist$be_best < 1.4), "| 1.4 to 2.8 (evidence range):",
    sum(dist$be_best >= 1.4 & dist$be_best < 2.8), "| 2.8 or more:", sum(dist$be_best >= 2.8), "\n")
cat("  realistic setting, below 1:", sum(dist$be_real < 1), "| zero (cannot reach 90% even with perfect timeliness):",
    sum(dist$be_real == 0), "\n")
print(dist[, c("region", "district", "mean_q", "be_best", "be_best_lo", "be_best_hi", "be_real", "p_best_below1")], row.names = FALSE)

# ---- break-even campaign coverage -----------------------------------------------------
# The ratio scales as 1 / C_c, so routine delivery reaches the 90% line whenever the campaign's
# per-cycle coverage is below C_c* = (routine coverage, burden-weighted) / 0.9. Only Upper East's
# campaign coverage is measured (87%); elsewhere 87% is assumed. Computed at k = 1.4 (lower end of
# the evidence) and k = 2 (central), best case and realistic (age gap + adherence), over the
# step-2 draws (and age-ratio draws for the realistic setting).
ccstar <- function(d_rows, k, setting, T = GREEN) sapply(seq_len(N_DRAWS), function(j) {
  ra <- if (setting == "best") 1 else ra_d[j]; ma <- if (setting == "best") 1 else madh_real
  r <- dreg[d_rows[1]]
  rk <- ratio_k(colSums(S[d_rows, , j, drop = FALSE]), colSums(Sd[d_rows, , j, drop = FALSE]),
                sum(A[d_rows]), sum(Cs[d_rows]), P1[r], 1, ra, ma, 1)
  interp_k(rk, k) / T
})
CCS <- list(k14_best = c(1.4, "best"), k2_best = c(2, "best"), k2_real = c(2, "real"))
for (nm in names(CCS)) {
  M <- sapply(seq_len(nD), function(d) ccstar(d, as.numeric(CCS[[nm]][1]), CCS[[nm]][2]))
  dist[[paste0("cc_", nm)]] <- round(apply(M, 2, median), 3)[match(dist$district, adm2$shapeName)]
  if (nm == "k2_best") {
    dist$cc_k2_best_lo <- round(apply(M, 2, quantile, 0.025), 3)[match(dist$district, adm2$shapeName)]
    dist$cc_k2_best_hi <- round(apply(M, 2, quantile, 0.975), 3)[match(dist$district, adm2$shapeName)]
  }
}
write.csv(dist, "data/derived/district_breakeven.csv", row.names = FALSE)
regcc <- do.call(rbind, lapply(NORTH, function(r) {
  i <- which(dreg == r)
  q <- function(k, s) quantile(ccstar(i, k, s), c(0.5, 0.025, 0.975))
  a <- q(1.4, "best"); b <- q(2, "best"); c2 <- q(2, "real")
  data.frame(region = r, cc_k14_best = a[1], cc_k2_best = b[1], cc_k2_best_lo = b[2], cc_k2_best_hi = b[3],
             cc_k2_real = c2[1], cc_k2_best_red = median(ccstar(i, 2, "best", RED)))
}))
regcc[, -1] <- round(regcc[, -1], 3)
write.csv(regcc, "data/derived/region_breakeven_cc.csv", row.names = FALSE)
cat("\nBreak-even campaign coverage per cycle (routine reaches 90% of the campaign if the campaign is below this):\n")
print(regcc, row.names = FALSE)
cat("Districts where the k = 2 best-case break-even coverage is below 0.75 (the low end assumed outside Upper East):",
    sum(dist$cc_k2_best < 0.75), "of", nD, "\n")

png("out/breakeven_cc.png", width = 1800, height = 950, res = 170)
par(mar = c(4.5, 8, 3.5, 1))
yy <- rev(seq_along(NORTH))
plot(NA, xlim = c(0.4, 1), ylim = c(0.2, length(NORTH) + 0.5), yaxt = "n", xlab = "Campaign coverage per cycle",
     ylab = "", main = "How low would campaign coverage have to fall for routine delivery to reach 90% of it?",
     cex.main = 0.95)
rect(CC_OTHER["lo"], 0, CC_OTHER["hi"], length(NORTH) + 1, col = "grey92", border = NA)
abline(v = CC_UE["mean"], lty = 2, col = "grey30")
text(CC_UE["mean"], length(NORTH) + 0.45, "87% measured (Upper East)", cex = 0.7, pos = 2, col = "grey30")
text(CC_OTHER["lo"], 0.3, "range assumed elsewhere", cex = 0.7, pos = 4, col = "grey40")
axis(2, yy, NORTH, las = 1, cex.axis = 0.85)
segments(regcc$cc_k2_best_lo, yy, regcc$cc_k2_best_hi, yy, lwd = 2, col = "#1b6ca8")
points(regcc$cc_k2_best, yy, pch = 16, col = "#1b6ca8", cex = 1.3)
points(regcc$cc_k14_best, yy, pch = 1, col = "#1b6ca8", cex = 1.3)
points(regcc$cc_k2_real, yy, pch = 17, col = "#c0392b", cex = 1.2)
legend("bottomleft", c("k = 2, best case (95% range)", "k = 1.4, best case", "k = 2, with age gap and adherence"),
       pch = c(16, 1, 17), col = c("#1b6ca8", "#1b6ca8", "#c0392b"), bty = "n", cex = 0.75)
dev.off()

# ---- favourable-to-routine check --------------------------------------------------
fav <- function(cc_fun, k = 1) {
  r <- sapply(seq_len(nD), function(d) mean(sapply(seq_len(N_DRAWS), function(j)
    interp_k(ratio_k(S[d, , j], Sd[d, , j], A[d], Cs[d], P1[dreg[d]], cc_fun(dreg[d]), 1, 1, 1), k))))
  c(green = sum(r >= GREEN), amber = sum(r >= RED & r < GREEN), red = sum(r < RED),
    pct_children_not_green = round(100 * sum(dist$children_3_59m[match(adm2$shapeName, dist$district)][r < GREEN]) / sum(dist$children_3_59m)))
}
fv <- rbind(`k=1, all else favourable, C_c mean` = fav(function(r) ccpar(r)["mean"]),
            `k=1, all else favourable, C_c low` = fav(function(r) ccpar(r)["lo"]),
            `k=1.4 (evidence lower bound), favourable, C_c mean` = fav(function(r) ccpar(r)["mean"], 1.4))
cat("\nFavourable-to-routine check:\n"); print(fv)
write.csv(data.frame(setting = rownames(fv), fv), "data/derived/favourable_check.csv", row.names = FALSE)

# ---- global sensitivity -------------------------------------------------------------
NG <- 20000
set.seed(20261022)
X <- data.frame(draw = rep(seq_len(N_DRAWS), length.out = NG),
                k = exp(rnorm(NG, log(K_MEDIAN), K_SDLOG)),
                cc_ue = rbeta_ci(NG, CC_UE["mean"], CC_UE["lo"], CC_UE["hi"]),
                cc_other = rbeta_ci(NG, CC_OTHER["mean"], CC_OTHER["lo"], CC_OTHER["hi"]),
                ra = runif(NG, AGE$lo, 1), a_r = runif(NG, ADH_R, ADH_C), f = runif(NG, 0, 1),
                reach = runif(NG, 0.9, 1))
RS  <- sapply(NORTH, function(r) colSums(S[dreg == r, , , drop = FALSE]), simplify = "array")   # k x draw x region
RSd <- sapply(NORTH, function(r) colSums(Sd[dreg == r, , , drop = FALSE]), simplify = "array")
RA <- tapply(A, factor(dreg, NORTH), sum); RC <- tapply(Cs, factor(dreg, NORTH), sum)
RQ <- sapply(NORTH, function(r) colSums(QB[dreg == r, , drop = FALSE]) / RA[r])                  # draw x region
saveRDS(list(KG = KG, RS = RS, RSd = RSd, RA = RA, RC = RC, RQ = RQ, P1 = P1, NORTH = NORTH),
        "data/derived/gsa_region_sums.rds")                     # for 33_sobol_total.R
kidx <- pmin(findInterval(X$k, KG), nK - 1); wk <- (X$k - KG[kidx]) / (KG[kidx + 1] - KG[kidx])
Y <- sapply(NORTH, function(r) {
  s1 <- RS[cbind(kidx, X$draw, match(r, NORTH))] * (1 - wk) + RS[cbind(kidx + 1, X$draw, match(r, NORTH))] * wk
  s2 <- RSd[cbind(kidx, X$draw, match(r, NORTH))] * (1 - wk) + RSd[cbind(kidx + 1, X$draw, match(r, NORTH))] * wk
  cc <- if (r == "Upper East") X$cc_ue else X$cc_other
  X$reach * P1[r] * madh_f(X$a_r, X$f) / cc * (s1 - (1 - X$ra) * s2) / RA[r]
})
first_order <- function(x, y, nb = 20) {
  b <- cut(x, unique(quantile(x, seq(0, 1, length.out = nb + 1))), include.lowest = TRUE)
  var(tapply(y, b, mean)[as.integer(b)]) / var(y)
}
FAC <- c(k = "k (routine vs vaccination slip)", cc = "campaign coverage", ra = "age reach (24-59 m)",
         a_r = "routine adherence", f = "partial-course protection", reach = "first-contact reach",
         q = "step-2 timeliness surface")
gsa <- do.call(rbind, lapply(NORTH, function(r) {
  y <- Y[, r]
  xs <- list(k = X$k, cc = if (r == "Upper East") X$cc_ue else X$cc_other, ra = X$ra, a_r = X$a_r,
             f = X$f, reach = X$reach, q = RQ[X$draw, r])
  data.frame(region = r, factor = FAC[names(xs)], first_order = round(sapply(xs, first_order, y = y), 3))
}))
write.csv(gsa, "data/derived/gsa_indices.csv", row.names = FALSE)
rr <- data.frame(region = NORTH, ratio_median = round(apply(Y, 2, median), 3),
                 lo = round(apply(Y, 2, quantile, 0.025), 3), hi = round(apply(Y, 2, quantile, 0.975), 3),
                 p_green = round(colMeans(Y >= GREEN), 3), p_red = round(colMeans(Y < RED), 3))
write.csv(rr, "data/derived/gsa_region_ratio.csv", row.names = FALSE)
cat("\nGlobal sensitivity: region ratio under all uncertainties together:\n"); print(rr, row.names = FALSE)
cat("\nFirst-order indices (share of variance in the region ratio):\n")
print(reshape(gsa, idvar = "factor", timevar = "region", direction = "wide"), row.names = FALSE)

# ---- figures ----------------------------------------------------------------------------
nr <- adm1[sub(" Region$", "", adm1$shapeName) %in% NORTH, ]
brk <- c(-Inf, 1, 1.4, 2.8, Inf)
lab <- c("< 1: routine would have to beat vaccination", "1 - 1.4: below the evidence", "1.4 - 2.8: within the evidence", "> 2.8: above it")
pal <- c("#8e1b1b", "#d9774f", "#f2cc7a", "#6cb36a")
png("out/breakeven_map.png", width = 2400, height = 1350, res = 170)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 3, 0.5), oma = c(4, 0, 4.5, 0))
for (v in c("be_best", "be_best_red")) {
  x <- pmin(dist[[v]][match(adm2$shapeName, dist$district)], 99)   # Inf: never crosses on the k grid
  plot(st_geometry(adm2), col = pal[cut(x, brk, right = FALSE)], border = "white", lwd = 0.5,
       main = c(be_best = "k at which routine falls below 90% of the campaign",
                be_best_red = "k at which routine falls below 70% of the campaign")[v], cex.main = 1)
  plot(st_geometry(nr), add = TRUE, lwd = 1.2)
}
mtext(paste0("Break-even k: how much more often a routine SMC visit can be missed than a vaccination visit\n",
             "Best case for routine delivery (no age gap, no adherence loss); district median. Evidence puts k near 2 (1.4-2.8)."),
      side = 3, outer = TRUE, line = 0.8, cex = 0.95)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new(); legend("bottom", lab, fill = pal, border = NA, bty = "n", horiz = TRUE, cex = 0.8, title = "break-even k")
dev.off()

png("out/gsa_chart.png", width = 1900, height = 1000, res = 170)
par(mar = c(4, 13, 3, 1))
M <- sapply(NORTH, function(r) gsa$first_order[gsa$region == r]); rownames(M) <- FAC
M <- M[order(rowMeans(M)), ]
barplot(t(M), beside = TRUE, horiz = TRUE, las = 1, xlim = c(0, 1), cex.names = 0.8,
        col = c("#1b6ca8", "#3a9e5f", "#c77c11", "#8e44ad", "#c0392b"),
        xlab = "Share of variance in the routine : campaign ratio (first-order)",
        main = "What drives the answer, by region", cex.main = 1)
legend("bottomright", rev(NORTH), fill = rev(c("#1b6ca8", "#3a9e5f", "#c77c11", "#8e44ad", "#c0392b")), bty = "n", cex = 0.8)
dev.off()
cat("done\n")
