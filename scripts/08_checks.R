# Checks on the step 2 timeliness model (06_inla_timeliness.R):
#   A. why population-weighted regional predictions sit below the observed shares
#   B. spatial block cross-validation against simpler models
#   C. sensitivity: region effects; motorised instead of walking time
# Writes data/derived/cv_results.csv and prints everything for the worklog.

source("scripts/lib_inla.R")
set.seed(20261001)
cl <- load_clusters()
mesh <- make_mesh(cl)
gd <- north_grid()
nth <- cl$region %in% NORTH
pct <- function(x) round(100 * x)
wmedian <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) >= sum(w) / 2)[1]] }

MAIN <- list(tt = "motor", spatial = TRUE)   # main model since 1 Oct 2026

# ---- A. decomposition ----------------------------------------------------------
m <- fit_model(cl, mesh, MAIN)
Pc <- draw_p(m, cl[nth, ])                       # new cluster at each surveyed location
Pg <- draw_p(m, gd)

A <- cl[nth, ] |> mutate(pc = rowMeans(Pc), fitted = m$p[nth]) |>
  group_by(region) |>
  summarise(observed = pct(sum(y) / sum(n)),
            observed_dhs_wt = pct(sum(weight * y) / sum(weight * n)),
            fitted = pct(sum(fitted * n) / sum(n)),
            pred_at_clusters = pct(sum(pc * n) / sum(n)),
            med_walk_children = round(median(rep(tt_walk_any, n))),
            med_motor_children = round(median(rep(tt_motor_any, n))))
G <- gd |> mutate(p = rowMeans(Pg)) |> group_by(region) |>
  summarise(pred_u5_weighted = pct(sum(p * w) / sum(w)),
            med_walk_u5 = round(wmedian(tw, w)), med_motor_u5 = round(wmedian(tm, w)))
cat("A. Northern regions, % of monthly contacts slipped (per contact)\n")
print(as.data.frame(inner_join(A, G, by = "region")), row.names = FALSE)

# ---- B. spatial block cross-validation -----------------------------------------
# Ten spatial blocks (k-means on cluster coordinates). Each block is held out in
# turn; its clusters are predicted as new clusters (fresh iid effect).
blocks <- kmeans(cbind(cl$X, cl$Y), centers = 10, nstart = 20)$cluster
MODELS <- list(
  intercept_only = list(tt = "motor", covars = FALSE, spatial = FALSE),
  covariates     = list(tt = "motor", spatial = FALSE),
  main_spatial   = MAIN,
  plus_region    = list(tt = "motor", spatial = TRUE, region = TRUE),
  walking        = list(tt = "walk", spatial = TRUE))
cv <- do.call(rbind, lapply(names(MODELS), function(nm) {
  p <- rep(NA_real_, nrow(cl))
  for (b in 1:10) {
    test <- which(blocks == b)
    p[test] <- fit_model(cl, mesh, MODELS[[nm]], test)$p[test]
  }
  y <- cl$y; n <- cl$n
  score <- function(k) c(
    logscore = -sum(dbinom(y[k], n[k], p[k], log = TRUE)) / sum(n[k]),
    brier = sum(y[k] * (1 - p[k])^2 + (n[k] - y[k]) * p[k]^2) / sum(n[k]),
    corr = cov.wt(cbind(y[k] / n[k], p[k]), wt = n[k], cor = TRUE)$cor[1, 2])
  data.frame(model = nm, subset = c("all Ghana", "north"),
             rbind(score(seq_along(y)), score(which(nth))))
}))
cv[, 3:5] <- round(cv[, 3:5], 4)
write.csv(cv, "data/derived/cv_results.csv", row.names = FALSE)
cat("\nB. Spatial block CV (10 blocks). Lower logscore/brier is better; corr = child-weighted\n",
    "   correlation of observed and predicted cluster shares.\n")
print(cv, row.names = FALSE)

# ---- C. sensitivity: odds ratio and regional predictions -----------------------
cat("\nC. Sensitivity (full data). OR per doubling of travel time; u5-weighted regional %\n")
for (nm in c("main_spatial", "plus_region", "walking")) {
  mm <- fit_model(cl, mesh, MODELS[[nm]])
  fx <- mm$fit$summary.fixed["x1", ]
  rs <- region_summary(draw_p(mm, gd, 200), gd$region, gd$w)
  cat(sprintf("  %-13s OR %.2f (%.2f-%.2f) | WAIC %.1f\n", nm, exp(fx$mean),
              exp(fx$`0.025quant`), exp(fx$`0.975quant`), mm$fit$waic$waic))
  cat("    ", paste(names(rs), rs, sep = " ", collapse = " | "), "\n")
}
