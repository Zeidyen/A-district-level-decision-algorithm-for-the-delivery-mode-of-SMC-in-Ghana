# District-level check of the step 2 spatial model against simpler models.
#
# Question a reviewer will ask: does the INLA-SPDE model predict district slip rates better than
# regional means? Step 4 uses the model at district level, so it is tested there.
#
# Leave-one-district-out: for each northern district, its DHS clusters are hidden, each model is
# refitted on all other clusters in Ghana, and the district's pooled slip rate is predicted.
# The spatial block CV in 08_checks.R holds out whole blocks, which can remove most of a region
# and so is unfair to region effects; holding out one district leaves the rest of its region in.
#
# Models (same code as step 2, lib_inla.R):
#   intercept   one national rate
#   region      regional means (region effects only; the "regiononly" step 4 variant)
#   covariates  travel time + population density, no spatial field (the "nospatial" variant)
#   main        covariates + Matern spatial field (the main model)
#   main_region main + region effects
#
# Scores, over held-out districts (contacts as weights):
#   logscore  -sum log Binomial(y | n, p) over held-out clusters / contacts (lower is better)
#   mae       contact-weighted mean |predicted - observed| district slip rate
#   corr      contact-weighted correlation of predicted and observed district slip rates
# plus paired comparisons with the regional-means model (districts where each model is closer,
# and a Wilcoxon signed-rank test on the district absolute errors).
#
# Clusters are assigned to districts by their displaced DHS coordinates (up to 2 km urban,
# 5 km rural), so a few sit in a neighbouring district; this affects all models alike.
#
# Outputs: data/derived/district_cv.csv (per district and model), district_cv_summary.csv,
#          out/district_cv.png

source("scripts/lib_inla.R")
suppressMessages(library(parallel))
cl <- load_clusters("slip56")
mesh <- make_mesh(cl)

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)])
adm2 <- adm2[adm2$region %in% NORTH, ]
cl$district <- adm2$shapeName[sapply(st_intersects(st_as_sf(cl, coords = c("LONGNUM", "LATNUM"), crs = 4326), adm2),
                                     function(i) if (length(i)) i[1] else NA)]
dists <- sort(unique(na.omit(cl$district)))
cat(sum(!is.na(cl$district)), "northern clusters in", length(dists), "districts; clusters per district:",
    paste(range(table(cl$district)), collapse = "-"), "\n")

MODELS <- list(
  intercept   = list(tt = "motor", spatial = FALSE, covars = FALSE, region = FALSE),
  region      = list(tt = "motor", spatial = FALSE, covars = FALSE, region = TRUE),
  covariates  = list(tt = "motor", spatial = FALSE, covars = TRUE,  region = FALSE),
  main        = list(tt = "motor", spatial = TRUE,  covars = TRUE,  region = FALSE),
  main_region = list(tt = "motor", spatial = TRUE,  covars = TRUE,  region = TRUE))

jobs <- expand.grid(model = names(MODELS), district = dists, stringsAsFactors = FALSE)
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(i) {
  test <- which(cl$district == jobs$district[i])
  m <- fit_model(cl, mesh, MODELS[[jobs$model[i]]], test = test)
  p <- m$p[test]
  data.frame(model = jobs$model[i], district = jobs$district[i], region = as.character(cl$region[test[1]]),
             clusters = length(test), n = sum(cl$n[test]), y = sum(cl$y[test]),
             pred = sum(cl$n[test] * p) / sum(cl$n[test]),
             ll = sum(dbinom(cl$y[test], cl$n[test], p, log = TRUE)))
}, mc.cores = 4)
bad <- !sapply(res, is.data.frame)
if (any(bad)) stop("failed fits: ", paste(jobs$model[bad], jobs$district[bad], collapse = "; "))
res <- do.call(rbind, res)
res$obs <- res$y / res$n
cat("fits done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
write.csv(res, "data/derived/district_cv.csv", row.names = FALSE)

wcor <- function(x, y, w) cov.wt(cbind(x, y), w / sum(w), cor = TRUE)$cor[1, 2]
summ <- do.call(rbind, lapply(names(MODELS), function(nm) {
  r <- res[res$model == nm, ]
  data.frame(model = nm, districts = nrow(r), logscore = -sum(r$ll) / sum(r$n),
             mae = sum(r$n * abs(r$pred - r$obs)) / sum(r$n), corr = wcor(r$pred, r$obs, r$n),
             corr_within_region = {
               e <- function(v) v - ave(v * r$n, r$region, FUN = sum) / ave(r$n, r$region, FUN = sum)
               wcor(e(r$pred), e(r$obs), r$n)
             })
}))
# Paired against regional means.
rg <- res[res$model == "region", ]
summ$closer_than_region <- sapply(summ$model, function(nm) {
  r <- res[res$model == nm, ]; r <- r[match(rg$district, r$district), ]
  sum(abs(r$pred - r$obs) < abs(rg$pred - rg$obs))
})
summ$wilcoxon_p_vs_region <- sapply(summ$model, function(nm) {
  if (nm == "region") return(NA)
  r <- res[res$model == nm, ]; r <- r[match(rg$district, r$district), ]
  wilcox.test(abs(r$pred - r$obs), abs(rg$pred - rg$obs), paired = TRUE)$p.value
})
summ[, 3:8] <- round(summ[, 3:8], 4)
write.csv(summ, "data/derived/district_cv_summary.csv", row.names = FALSE)
cat("\nLeave-one-district-out CV, northern districts (lower logscore / mae is better):\n")
print(summ, row.names = FALSE)

# Noise floor: with this few contacts per district, how much would the observed rate scatter
# around a perfectly known true rate? Expected |error| under binomial sampling at the main prediction.
mn <- res[res$model == "main", ]
floor_mae <- sum(mn$n * sqrt(2 / pi) * sqrt(mn$pred * (1 - mn$pred) / mn$n)) / sum(mn$n)
cat(sprintf("\nBinomial noise floor (expected MAE if predictions were exactly right, no overdispersion): %.4f\n", floor_mae))
cat("median contacts per district:", median(mn$n), "| clusters:", median(mn$clusters), "\n")
write.csv(data.frame(noise_floor_mae = round(floor_mae, 4), median_contacts = median(mn$n)),
          "data/derived/district_cv_noise.csv", row.names = FALSE)

# Figure: predicted vs observed, regional means vs main model.
png("out/district_cv.png", width = 2000, height = 1000, res = 170)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
cols <- setNames(c("#1b6ca8", "#3a9e5f", "#c77c11", "#8e44ad", "#c0392b"), NORTH)
for (nm in c("region", "main")) {
  r <- res[res$model == nm, ]; s <- summ[summ$model == nm, ]
  plot(r$pred, r$obs, cex = sqrt(r$n) / 6, pch = 21, bg = adjustcolor(cols[as.character(r$region)], 0.6), col = "grey30",
       xlim = c(0, 0.5), ylim = c(0, 0.5), xlab = "Predicted district slip rate (held out)",
       ylab = "Observed district slip rate (DHS)",
       main = sprintf("%s: MAE %.3f, corr %.2f", c(region = "Regional means", main = "Main spatial model")[nm],
                      s$mae, s$corr), cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey50")
  if (nm == "region") legend("topleft", NORTH, pt.bg = cols, pch = 21, bty = "n", cex = 0.75)
}
dev.off()
cat("done\n")
