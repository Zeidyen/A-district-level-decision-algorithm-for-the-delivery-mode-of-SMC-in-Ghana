# SMC delivery-mode algorithm: a district-level decision rule for campaign vs routine delivery of
# seasonal malaria chemoprevention, designed to sit inside a national malaria programme's annual
# planning cycle and to be re-run each year as new data arrive.
#
# Per SMC-eligible district (eligibility is decided upstream by the national SMC criteria):
#   1 CAMPAIGN CHECK   where campaign coverage per cycle is measured BY HOUSEHOLD SURVEY (end-of-cycle
#                      or end-of-season coverage survey): is it below the programme's target (CC_TARGET)?
#                      yes -> "Fix the campaign first": a weak campaign makes any comparison misleading.
#                      Not measured -> flagged, 87% assumed. Administrative coverage (doses distributed /
#                      target population) is usually 95-100% or more and cannot tell a good campaign from
#                      a weak one, so it is never used here as measured coverage.
#   2 AGE REACH        which routine model is proposed? "clinic" (EPI-linked contacts, which barely
#                      reach children aged 2-4: the age gap and adherence loss apply) or "outreach"
#                      (built to reach 2-4-year-olds; assumed, to be verified by the pilot).
#   3 RELIABILITY      break-even k (how much more often a routine SMC visit can be missed than a
#                      vaccination visit before routine delivery falls below 90% of the campaign),
#                      with its uncertainty, against the evidence range for k (K_LO-K_HI):
#                        P(break-even < K_LO) >= P_CONF          -> "Keep campaign"
#                        P(break-even > K_HI) >= P_CONF          -> switch candidate (step 4)
#                        median < K_LO but not confidently        -> "Measure first"
#                        otherwise (median in or above the range) -> "Pilot"
#   4 EQUITY           for Pilot and switch candidates: if the poorest wealth quintile bears a
#                      significantly larger share of the protection lost than its share of children
#                      (disproportionality ratio = loss share / child share; lower 95% bootstrap bound
#                      above EQUITY_REF = 1, i.e. losses not proportional to population; region level,
#                      as districts have too few card-holding children), a switch candidate becomes
#                      "Hybrid" (routine plus a campaign mop-up for the poorest and most remote), else
#                      "Switch"; pilots are flagged to test a safeguard.
#   5 PRIORITY         pilot districts ranked by routine-to-campaign ratio (the hybrid curve order)
#   ADMIN RESOLUTION   a "Measure first" district probably keeps the campaign but rests on unmeasured
#                      campaign coverage. If administrative coverage is available AND can be corrected
#                      against a survey (ADMIN_CALIB = survey / administrative coverage where both exist,
#                      e.g. Upper East), the corrected figure resolves it: corrected coverage at or above
#                      the district's break-even campaign coverage at k = K_LO (the coverage below which
#                      routine delivery could match the campaign even at the low end of the evidence for k)
#                      -> "Keep campaign". Inflated administrative coverage can only make routine
#                      delivery look worse, never better, so uncorrected figures are not used.
#   6 MEASURE          for Pilot, Measure first, Switch and Hybrid districts: what to measure first,
#                      ordered by the global sensitivity analysis for the district's region.
#                      After a pilot (rule, not computed here): switch if the measured k is below the
#                      district's break-even and age reach and adherence are confirmed; otherwise revert.
# Robustness flag: districts whose outcome disagrees with their class across the 14 step 4 variants
# (a pilot or switch district never green in any variant; a keep-campaign district green in some).
#
# Inputs: district_breakeven.csv and gsa_indices.csv (21), district_step4.csv (10),
#         district_robustness.csv (14), equity_by_group.csv (30).
# Outputs: data/derived/delivery_algorithm.csv, delivery_algorithm_counts.csv;
#          out/delivery_algorithm_map.png, out/delivery_algorithm_flow.png

suppressMessages({library(sf); library(dplyr)})
sf_use_s2(FALSE)

# ---- settings (programme choices; change here) -----------------------------------------
K_LO <- 1.4; K_HI <- 2.8                       # 95% range of the evidence for k (lognormal, median 2)
P_CONF <- 0.8                                  # probability needed to keep the campaign or switch
EQUITY_REF <- 1                                # safeguard if the ratio's lower 95% bound exceeds this
# (no WHO or literature threshold exists for this ratio; 1 means losses proportional to population)
MEASURED_CC <- c(`Upper East` = 0.87)          # measured campaign coverage per cycle (Ayamba 2025)
# Three campaign-coverage quantities, with different roles: CC_TARGET is a quality gate (step 1) on
# survey-measured coverage; C_c, the coverage in the routine-to-campaign ratio and break-even k (step 3),
# is the survey value where measured and 87% (75-95%) otherwise -- set in 21_breakeven_gsa.R (CC_UE,
# CC_OTHER), which must agree with MEASURED_CC; break-even campaign coverage is an output.
CC_TARGET <- 0.80                              # programme target for campaign coverage per cycle
# Administrative SMC coverage from NMEP, if supplied: columns region, district, coverage (per cycle,
# averaged over cycles and recent years). Calibration: survey / administrative coverage in regions with
# both; NA (no correction possible) means administrative figures are not used.
ADMIN_FILE <- "data/raw/nmep/smc_admin_coverage.csv"
admin <- if (file.exists(ADMIN_FILE)) read.csv(ADMIN_FILE) else NULL
ADMIN_CALIB <- if (is.null(admin)) NA else {
  both <- intersect(names(MEASURED_CC), admin$region)
  if (length(both)) mean(MEASURED_CC[both] / tapply(admin$coverage, admin$region, mean)[both]) else NA
}
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
clean <- function(x) trimws(gsub("\\s+", " ", x))

be <- read.csv("data/derived/district_breakeven.csv")
s4 <- read.csv("data/derived/district_step4.csv")
rb <- read.csv("data/derived/district_robustness.csv")
eq <- read.csv("data/derived/equity_region_ratio_ci.csv") |>
  mutate(region = tools::toTitleCase(region), eq_ratio = ratio, eq_lo = lo)
gsa <- read.csv("data/derived/gsa_indices.csv")
be$district <- clean(be$district); s4$district <- clean(s4$district); rb$district <- clean(rb$district)
be$ratio <- s4$ratio[match(be$district, s4$district)]
be$eq_ratio <- eq$eq_ratio[match(be$region, eq$region)]
be$eq_lo <- eq$eq_lo[match(be$region, eq$region)]
rb$grp <- with(rb, ifelse(red == n, "red in all variants", ifelse(red >= 9, "red in most variants",
                     ifelse(green == 0, "never green", "green in some variants"))))
be$robust_grp <- rb$grp[match(be$district, rb$district)]

# What to measure first, by region: inputs ordered by their share of variance in the region's
# routine-to-campaign ratio (the timeliness surface is part of routine contact reliability).
MEAS <- c("k (routine vs vaccination slip)" = "routine contact reliability (k)",
          "step-2 timeliness surface" = "routine contact reliability (k)",
          "campaign coverage" = "campaign coverage", "age reach (24-59 m)" = "reach of children aged 2-4",
          "first-contact reach" = "first-contact reach", "routine adherence" = "routine adherence",
          "partial-course protection" = "partial-course protection")
measure <- sapply(NORTH, function(r) {
  g <- gsa[gsa$region == r, ]; g$m <- MEAS[g$factor]
  v <- tapply(g$first_order, g$m, sum); paste(names(sort(v, decreasing = TRUE))[1:3], collapse = "; ")
})

run <- function(model) {
  real <- model == "clinic"
  bek <- if (real) be$be_real else be$be_best
  plo <- if (real) be$p_real_below_klo else be$p_best_below_klo
  phi <- if (real) be$p_real_above_khi else be$p_best_above_khi
  meas <- MEASURED_CC[be$region]
  step1 <- ifelse(is.na(meas), "not measured (87% assumed)", ifelse(meas < CC_TARGET, "below target", "meets target"))
  # Steps in the order of the flow diagram: 1 campaign check; 2 age reach (a clinic/EPI-linked model
  # does not reach children aged 2-4 -> keep campaign by rule; the break-even with the age gap applied
  # is kept as a quantitative confirmation, `step2_confirmed`); 3 reliability.
  step2 <- rep(if (real) "does not reach 2-4" else "to be verified", nrow(be))
  out <- ifelse(step1 == "below target", "Fix the campaign first",
         ifelse(step2 == "does not reach 2-4", "Keep campaign",
         ifelse(plo >= P_CONF, "Keep campaign", ifelse(phi >= P_CONF, "Switch",
         ifelse(bek < K_LO, "Measure first", "Pilot")))))
  step2_confirmed <- if (real) plo >= P_CONF else rep(NA, nrow(be))
  eqflag <- out %in% c("Pilot", "Switch") & be$eq_lo > EQUITY_REF
  out[out == "Switch" & eqflag] <- "Hybrid"
  # Administrative resolution of "Measure first" (needs corrected administrative coverage).
  adm <- if (is.null(admin) || is.na(ADMIN_CALIB)) rep(NA, nrow(be)) else
    ADMIN_CALIB * admin$coverage[match(be$district, clean(admin$district))]
  thr <- if (real) be$cc_k2_real else be$cc_k14_best
  resolved <- out == "Measure first" & !is.na(adm) & adm >= thr
  out[resolved] <- "Keep campaign"
  robust <- ifelse(out %in% c("Pilot", "Switch", "Hybrid") & be$robust_grp != "green in some variants",
                   paste("check:", be$robust_grp), ifelse(out == "Keep campaign" & be$robust_grp == "green in some variants",
                   "check: green in some variants", ""))
  d <- data.frame(model, region = be$region, district = be$district, children = be$children_3_59m,
                  breakeven_k = bek, p_below_klo = plo, p_above_khi = phi, campaign_check = step1,
                  age_reach = step2, step2_confirmed = step2_confirmed,
                  equity_ratio = round(be$eq_ratio, 2), equity_ratio_lo = round(be$eq_lo, 2), equity_flag = eqflag, outcome = out, ratio = be$ratio,
                  robustness_flag = robust, admin_cov_corrected = round(adm, 3), cc_threshold = round(thr, 3),
                  resolved_by_admin = resolved,
                  measure_first = ifelse(out %in% c("Pilot", "Measure first", "Switch", "Hybrid"), measure[be$region], ""))
  d$pilot_priority <- NA
  p <- which(d$outcome == "Pilot"); d$pilot_priority[p] <- rank(-d$ratio[p], ties.method = "first")
  d
}
OUT <- c("Keep campaign", "Measure first", "Pilot", "Hybrid", "Switch", "Fix the campaign first")
res <- rbind(run("clinic"), run("outreach"))
res <- res[order(res$model, factor(res$outcome, rev(OUT)), res$pilot_priority, -res$ratio), ]
write.csv(res, "data/derived/delivery_algorithm.csv", row.names = FALSE)

cnt <- res |> group_by(model, outcome = factor(outcome, OUT), .drop = FALSE) |>
  summarise(districts = n(), children = sum(children), .groups = "drop") |>
  group_by(model) |> mutate(pct_children = round(100 * children / sum(children), 1)) |> ungroup()
write.csv(cnt, "data/derived/delivery_algorithm_counts.csv", row.names = FALSE)
cat("Algorithm outcomes (districts, % of children):\n"); print(as.data.frame(cnt), row.names = FALSE)
o <- res[res$model == "outreach", ]
cat("\nOutreach model, pilot districts in priority order:", paste(o$district[o$outcome == "Pilot"], collapse = ", "), "\n")
cat("Measure first:", paste(o$district[o$outcome == "Measure first"], collapse = ", "), "\n")
cat("Switch / Hybrid:", paste(o$district[o$outcome %in% c("Switch", "Hybrid")], collapse = ", "), "\n")
cat("Robustness flags:", paste(o$district[o$robustness_flag != ""], o$robustness_flag[o$robustness_flag != ""], sep = " - ", collapse = "; "), "\n")
cat("Administrative coverage:", if (is.null(admin)) "not supplied" else sprintf("supplied; calibration factor %s", round(ADMIN_CALIB, 3)),
    "| measure-first districts resolved by it:", sum(res$resolved_by_admin), "\n")
cat("Equity ratio by region:", paste(eq$region, round(eq$eq_ratio, 2), collapse = " | "), "\n")
cat("Measure first, by region:\n"); for (r in NORTH) cat(" ", r, ":", measure[r], "\n")

# ---- map: outcome by district, both routine models ---------------------------------------
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)])
adm2 <- adm2[adm2$region %in% NORTH, ]; nr <- adm1[sub(" Region$", "", adm1$shapeName) %in% NORTH, ]
key <- clean(adm2$shapeName)
COL <- c(`Keep campaign` = "#b03a2e", `Measure first` = "#c39bd3", Pilot = "#f0b44c", Hybrid = "#7fb3d5",
         Switch = "#2e8b57", `Fix the campaign first` = "#7f8c8d")
png("out/delivery_algorithm_map.png", width = 2400, height = 1300, res = 170)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 3.5, 0.5), oma = c(3.5, 0, 0, 0))
for (mdl in c("clinic", "outreach")) {
  r <- res[res$model == mdl, ]
  oc <- r$outcome[match(key, r$district)]
  plot(st_geometry(adm2), col = COL[oc], border = "white", lwd = 0.5,
       main = c(clinic = "Routine model linked to clinic/EPI contacts\n(does not reach children aged 2-4)",
                outreach = "Routine model built to reach children aged 2-4\n(outreach or community volunteers; to be verified)")[mdl],
       cex.main = 0.95)
  plot(st_geometry(nr), add = TRUE, lwd = 1.2)
  fl <- adm2[key %in% r$district[r$robustness_flag != ""], ]
  if (nrow(fl)) plot(st_geometry(fl), add = TRUE, border = "black", lwd = 1.6, lty = 2)
  pp <- r[r$outcome == "Pilot" & r$pilot_priority <= 6, ]
  if (nrow(pp)) {
    top <- adm2[match(pp$district, key), ]
    text(st_coordinates(st_point_on_surface(top)), labels = pp$pilot_priority, cex = 0.7, font = 2)
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
shown <- OUT[OUT %in% res$outcome]
plot.new(); legend("bottom", shown, fill = COL[shown], border = NA, bty = "n", horiz = TRUE, cex = 0.8,
                   title = "Outcome (numbers: pilot priority 1-6; dashed outline: robustness flag)")
dev.off()

# ---- flowchart with district counts (outreach model; clinic model in its own box) --------
# Black and white: decision steps white, outcomes light grey with bold text; large text for print.
n_o <- function(x) { z <- cnt[cnt$model == "outreach" & cnt$outcome == x, ]
  sprintf("%d district%s, %.0f%% of children", z$districts, ifelse(z$districts == 1, "", "s"), z$pct_children) }
n_c <- function(x) { z <- cnt[cnt$model == "clinic" & cnt$outcome == x, ]; z$districts }
INK <- "black"; OUTF <- "#eeeeee"
png("out/delivery_algorithm_flow.png", width = 2400, height = 3300, res = 220)
par(mar = c(0, 0, 0, 0)); plot(NA, xlim = c(-0.9, 11.0), ylim = c(-0.1, 15.0), axes = FALSE, xlab = "", ylab = "")
box <- function(x, y, w, h, lab, fill = "white", cex = 0.95, font = 1, lwd = 1.3) {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = INK, lwd = lwd)
  text(x, y, lab, cex = cex, font = font, col = INK)
}
LX <- 3.2; W <- 5.6; RX <- 9.15; RW <- 3.3
harr <- function(y, lab) {
  arrows(LX + W / 2, y, RX - RW / 2, y, length = 0.09, lwd = 1.3, col = INK)
  text((LX + W / 2 + RX - RW / 2) / 2, y + 0.05, lab, cex = 0.8, pos = 3, col = INK)
}
varr <- function(y0, y1, lab = NULL) {
  arrows(LX, y0, LX, y1, length = 0.09, lwd = 1.3, col = INK)
  if (!is.null(lab)) text(LX + 0.1, (y0 + y1) / 2, lab, cex = 0.8, pos = 4, col = INK)
}
H <- c(start = 0.9, s1 = 1.3, s2 = 1.3, s3 = 2.6, s4 = 2.0, s5 = 1.2, s6 = 1.5)
Y <- c(start = 14.3, s1 = 12.75, s2 = 10.9, s3 = 8.3, s4 = 5.35, s5 = 3.15, s6 = 1.2)
box(LX, Y[["start"]], W, H[["start"]], "SMC-eligible district\n(meets the national SMC eligibility criteria)", font = 2)
box(LX, Y[["s1"]], W, H[["s1"]], sprintf("1  Campaign check\nIs survey-measured campaign coverage\nper cycle below the programme target (%.0f%%)?", 100 * CC_TARGET))
box(LX, Y[["s2"]], W, H[["s2"]], "2  Age reach\nDoes the routine model reach children aged 2-4?\n(clinic/EPI contacts: no; outreach: to be verified)")
box(LX, Y[["s3"]], W, H[["s3"]], sprintf("3  Reliability\nBreak-even k: how much more often a routine\nSMC visit could be missed than a vaccination\nvisit before routine delivery falls below 90%%\nof the campaign. Compared, with its\nuncertainty, with the evidence for k (%.1f-%.1f)", K_LO, K_HI))
box(LX, Y[["s4"]], W, H[["s4"]], "4  Equity safeguard\nWould the poorest wealth quintile bear\nsignificantly more than its share of the\nprotection lost? (lower 95% bound of\nloss share / child share > 1)")
box(LX, Y[["s5"]], W, H[["s5"]], "5  Pilot priority\nRank pilot districts by routine-to-campaign\nratio; each pilot includes an equity safeguard")
box(LX, Y[["s6"]], W, H[["s6"]], "6  Measure and decide\nMeasure the region's main drivers first.\nSwitch if measured k is below break-even and\nage reach and adherence are confirmed; else revert")
yk <- Y[["s3"]] + 0.7; ym <- Y[["s3"]] - 0.7; yh <- Y[["s4"]] + 0.5; ysw <- Y[["s4"]] - 0.5
ob <- function(y, h, lab, cex = 0.9) box(RX, y, RW, h, lab, fill = OUTF, cex = cex, font = 2, lwd = 1.6)
ob(Y[["s1"]], 1.2, sprintf("Fix the campaign first\n%d districts\n(only Upper East measured:\n%.0f%%, above the target)", cnt$districts[cnt$model == "outreach" & cnt$outcome == "Fix the campaign first"], 100 * MEASURED_CC[["Upper East"]]), cex = 0.85)
ob(Y[["s2"]], 1.0, sprintf("Keep campaign\nclinic model:\n%d of 55 districts", n_c("Keep campaign")))
ob(yk, 1.0, paste0("Keep campaign\n", sub(", ", ",\n", n_o("Keep campaign"))))
ob(ym, 1.15, paste0("Measure first\n", n_o("Measure first"), "\n(survey or corrected\nadministrative coverage)"), cex = 0.82)
ob(yh, 0.9, paste0("Hybrid: routine +\ncampaign mop-up\n", n_o("Hybrid")), cex = 0.82)
ob(ysw, 0.9, paste0("Switch to routine\n", n_o("Switch")), cex = 0.85)
ob(Y[["s5"]], 1.0, paste0("Pilot\n", sub(", ", ",\n", n_o("Pilot"))))
varr(Y[["start"]] - H[["start"]] / 2, Y[["s1"]] + H[["s1"]] / 2)
varr(Y[["s1"]] - H[["s1"]] / 2, Y[["s2"]] + H[["s2"]] / 2, "no, or not measured (87% assumed in step 3)")
varr(Y[["s2"]] - H[["s2"]] / 2, Y[["s3"]] + H[["s3"]] / 2, "yes, or to be tested")
varr(Y[["s3"]] - H[["s3"]] / 2, Y[["s4"]] + H[["s4"]] / 2, sprintf("P(k* > %.1f) ≥ %.0f%%: switch candidate", K_HI, 100 * P_CONF))
varr(Y[["s5"]] - H[["s5"]] / 2, Y[["s6"]] + H[["s6"]] / 2)
harr(Y[["s1"]], "yes"); harr(Y[["s2"]], "no")
harr(yk, sprintf("P(k* < %.1f) ≥ %.0f%%", K_LO, 100 * P_CONF)); harr(ym, sprintf("k* probably < %.1f", K_LO))
harr(yh, "yes"); harr(ysw, "no"); harr(Y[["s5"]], "")
XL <- LX - W / 2 - 0.45
segments(LX - W / 2, Y[["s3"]] - 0.85, XL, Y[["s3"]] - 0.85, lwd = 1.3, col = INK)
segments(XL, Y[["s3"]] - 0.85, XL, Y[["s5"]], lwd = 1.3, col = INK)
arrows(XL, Y[["s5"]], LX - W / 2, Y[["s5"]], length = 0.09, lwd = 1.3, col = INK)
text(XL - 0.05, (Y[["s3"]] - 0.85 + Y[["s5"]]) / 2, sprintf("otherwise (k* likely within %.1f-%.1f)", K_LO, K_HI), cex = 0.8, srt = 90, col = INK, pos = 3)
text(5.05, 0.05, "k* = break-even k. Counts: northern Ghana, routine model built to reach children aged 2-4.", cex = 0.72, col = INK)
dev.off()
cat("done\n")
