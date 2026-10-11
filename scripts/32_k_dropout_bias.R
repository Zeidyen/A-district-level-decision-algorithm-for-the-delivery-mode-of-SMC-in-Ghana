# How biased is k when it is calibrated on dropout but applied to lateness?
#
# k scales the per-contact slip probability q, where a slip is a contact given more than 8 weeks
# after the previous dose OR never given. The evidence for k (median 2, 95% range 1.4-2.8) comes
# from malaria-vaccine vs pentavalent DROPOUT (doses never given). Split q into its parts:
#   q = q_never + q_late
# If a routine malaria-product contact raises both parts by the same factor, the effective k is
# the dropout ratio itself (the main analysis: upper bound). If it raised only the never-given part
# and were no later than a vaccination contact, the effective k would be
#   k_eff = (k q_never + q_late) / q = 1 + (k - 1) d,   d = q_never / q
# (lower bound). The truth lies between, and the main analysis is the conservative end: it can
# only overstate how much worse routine SMC is, not understate it. Districts whose break-even k is
# below 1 keep the campaign under any k >= 1, so they are unaffected.
#
# d is measured from card-dated Penta contacts in the north (27a_slip_features.R). Algorithm
# outcomes are recomputed with the evidence range mapped to the lower bound, approximating each
# district's break-even distribution as lognormal from its median and 95% interval.
#
# Outputs: data/derived/k_dropout_bias.csv (bounds), k_dropout_bias_outcomes.csv (district
#          outcomes under both bounds)

suppressMessages(library(dplyr))
d0 <- read.csv("data/derived/slip_contacts_features.csv") |> filter(north)
sl <- d0[d0$slip == 1, ]
d_all <- mean(sl$slip_type == "never given")
d_reg <- tapply(sl$slip_type == "never given", sl$region, mean)
keff <- function(k, d) 1 + (k - 1) * d
K <- c(lo = 1.4, med = 2, hi = 2.8)
bounds <- data.frame(scenario = c("upper bound (main: both parts scale)", "lower bound (only never-given part scales)"),
                     k_lo = c(K[["lo"]], keff(K[["lo"]], d_all)), k_med = c(K[["med"]], keff(K[["med"]], d_all)),
                     k_hi = c(K[["hi"]], keff(K[["hi"]], d_all)))
bounds[, 2:4] <- round(bounds[, 2:4], 2)
write.csv(rbind(bounds, data.frame(scenario = "never-given share of slips (north)", k_lo = round(d_all, 3), k_med = NA, k_hi = NA)),
          "data/derived/k_dropout_bias.csv", row.names = FALSE)
cat(sprintf("Never-given share of slipped contacts, north: %.3f (%s)\n", d_all,
            paste(names(d_reg), round(d_reg, 2), collapse = ", ")))
print(bounds, row.names = FALSE)

# Re-classify (outreach model) with the evidence range at each bound.
be <- read.csv("data/derived/district_breakeven.csv")
alg <- read.csv("data/derived/delivery_algorithm.csv") |> filter(model == "outreach")
be$district <- trimws(gsub("\\s+", " ", be$district))
mu <- log(pmax(be$be_best, 1e-3)); sdl <- (log(pmax(be$be_best_hi, 1e-3)) - log(pmax(be$be_best_lo, 1e-3))) / 3.92
cls <- function(klo, khi) {
  plo <- pnorm((log(klo) - mu) / pmax(sdl, 1e-6)); phi <- 1 - pnorm((log(khi) - mu) / pmax(sdl, 1e-6))
  ifelse(plo >= 0.8, "Keep campaign", ifelse(phi >= 0.8, "Switch candidate", ifelse(be$be_best < klo, "Measure first", "Pilot")))
}
out <- data.frame(region = be$region, district = be$district, breakeven_k = be$be_best,
                  upper = cls(bounds$k_lo[1], bounds$k_hi[1]), lower = cls(bounds$k_lo[2], bounds$k_hi[2]),
                  main_algorithm = alg$outcome[match(be$district, alg$district)])
write.csv(out, "data/derived/k_dropout_bias_outcomes.csv", row.names = FALSE)
cat("\nOutcomes (lognormal approximation; 'upper' reproduces the main algorithm):\n")
print(table(upper = out$upper, lower = out$lower))
cat("agreement of 'upper' with the main algorithm:", mean(out$upper == out$main_algorithm), "\n")
cat("districts with break-even k < 1 (unaffected by the bias):", sum(out$breakeven_k < 1), "\n")
cat("still 'Keep campaign' at the lower bound:", sum(out$lower == "Keep campaign"),
    "| children share:", round(sum(read.csv("data/derived/district_breakeven.csv")$children_3_59m[out$lower == "Keep campaign"]) /
                                sum(read.csv("data/derived/district_breakeven.csv")$children_3_59m), 3), "\n")
