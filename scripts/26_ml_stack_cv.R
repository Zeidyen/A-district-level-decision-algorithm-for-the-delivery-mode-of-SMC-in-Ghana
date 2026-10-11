# Does a stacked machine-learning covariate improve the step 2 spatial model?
#
# 25_ml_stack.py gives each DHS cluster a stacked prediction (mean logit of random forest,
# gradient boosting and ridge logistic on 16 satellite and access covariates), always from learners
# that did not see that cluster. Here it enters the INLA-SPDE model as a covariate:
#   stacked      ML covariate + Matern spatial field + cluster effect (Bhatt et al. 2017 style)
#   stacked_cov  as stacked, plus the main model's travel time and population density
#   ml_glm       ML covariate + cluster effect, no spatial field
# scored by leave-one-district-out over the northern districts, exactly as 22_district_cv.R, and
# compared with the models there and with the learners alone (no spatial field).
# Full-data fit: the stacked model on all clusters (10-fold out-of-fold ML covariate), with WAIC
# against the main model.
#
# Outputs: data/derived/ml_stack_cv.csv (per district and model), ml_stack_cv_summary.csv
#          (all models side by side), ml_stack_fixed.csv; out/ml_stack_cv.png

source("scripts/lib_inla.R")
suppressMessages(library(parallel))
cl <- load_clusters("slip56")
mesh <- make_mesh(cl)
full <- read.csv("data/derived/ml_oof_full.csv")
lodo <- read.csv("data/derived/ml_lodo.csv")
cl$ml <- full$ml_stack[match(cl$cluster, full$cluster)]
miss <- is.na(cl$ml)
cat("clusters without an ML prediction (missing covariates):", sum(miss), "\n")
cl <- cl[!miss, ]

# ---- full-data fits -------------------------------------------------------------------
SPEC <- list(main        = list(tt = "motor", spatial = TRUE,  covars = TRUE,  region = FALSE),
             stacked     = list(tt = "motor", spatial = TRUE,  covars = FALSE, region = FALSE, ml = TRUE),
             stacked_cov = list(tt = "motor", spatial = TRUE,  covars = TRUE,  region = FALSE, ml = TRUE),
             ml_glm      = list(tt = "motor", spatial = FALSE, covars = FALSE, region = FALSE, ml = TRUE))
fits <- lapply(SPEC, function(s) fit_model(cl, mesh, s))
fx <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]$fit; s <- f$summary.fixed
  data.frame(model = nm, term = rownames(s), mean = s$mean, lo = s$`0.025quant`, hi = s$`0.975quant`,
             waic = f$waic$waic)
}))
fx[, 3:6] <- round(fx[, 3:6], 3)
write.csv(fx, "data/derived/ml_stack_fixed.csv", row.names = FALSE)
cat("\nFull-data fits (WAIC lower is better; x3 = stacked ML covariate, logit scale):\n")
print(fx, row.names = FALSE)

# ---- leave-one-district-out ------------------------------------------------------------
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)])
adm2 <- adm2[adm2$region %in% NORTH, ]
cl$district <- adm2$shapeName[sapply(st_intersects(st_as_sf(cl, coords = c("LONGNUM", "LATNUM"), crs = 4326), adm2),
                                     function(i) if (length(i)) i[1] else NA)]
dists <- sort(unique(na.omit(cl$district)))
CVM <- c("stacked", "stacked_cov", "ml_glm")
jobs <- expand.grid(model = CVM, district = dists, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(i) {
  dd <- jobs$district[i]
  l <- lodo[lodo$held_out == dd, ]
  c2 <- cl; c2$ml <- l$ml[match(c2$cluster, l$cluster)]
  c2 <- c2[!is.na(c2$ml), ]
  test <- which(c2$district %in% dd)
  m <- fit_model(c2, mesh, SPEC[[jobs$model[i]]], test = test)
  p <- m$p[test]
  data.frame(model = jobs$model[i], district = dd, region = as.character(c2$region[test[1]]),
             clusters = length(test), n = sum(c2$n[test]), y = sum(c2$y[test]),
             pred = sum(c2$n[test] * p) / sum(c2$n[test]),
             ll = sum(dbinom(c2$y[test], c2$n[test], p, log = TRUE)))
}, mc.cores = 4)
bad <- !sapply(res, is.data.frame)
if (any(bad)) stop("failed fits: ", paste(jobs$model[bad], jobs$district[bad], collapse = "; "))
res <- do.call(rbind, res); res$obs <- res$y / res$n
write.csv(res, "data/derived/ml_stack_cv.csv", row.names = FALSE)

# All models side by side: spatial models (22_district_cv.R), learners alone (25), stacked (here).
old <- read.csv("data/derived/district_cv.csv")
lr <- read.csv("data/derived/ml_learner_cv.csv")
lr <- data.frame(model = paste0("ml_alone_", lr$learner), district = lr$district, n = lr$n, y = lr$y,
                 pred = lr$pred, ll = lr$ll, obs = lr$obs)
cols <- c("model", "district", "n", "y", "pred", "ll", "obs")
allr <- rbind(old[, cols], lr[, cols], res[, cols])
allr <- allr[allr$district %in% dists, ]
wcor <- function(x, y, w) cov.wt(cbind(x, y), w / sum(w), cor = TRUE)$cor[1, 2]
rg <- allr[allr$model == "region", ]; mn <- allr[allr$model == "main", ]
summ <- do.call(rbind, lapply(unique(allr$model), function(nm) {
  r <- allr[allr$model == nm, ]; a <- r[match(mn$district, r$district), ]
  data.frame(model = nm, districts = nrow(r), logscore = -sum(r$ll) / sum(r$n),
             mae = sum(r$n * abs(r$pred - r$obs)) / sum(r$n), corr = wcor(r$pred, r$obs, r$n),
             closer_than_main = sum(abs(a$pred - a$obs) < abs(mn$pred - mn$obs), na.rm = TRUE),
             wilcoxon_p_vs_main = if (nm == "main") NA else
               wilcox.test(abs(a$pred - a$obs), abs(mn$pred - mn$obs), paired = TRUE)$p.value)
}))
summ[, 3:7] <- round(summ[, 3:7], 4)
summ <- summ[order(summ$mae), ]
write.csv(summ, "data/derived/ml_stack_cv_summary.csv", row.names = FALSE)
cat("\nLeave-one-district-out, all models (lower log score / MAE is better):\n")
print(summ, row.names = FALSE)

png("out/ml_stack_cv.png", width = 1700, height = 950, res = 170)
par(mar = c(4.5, 12, 3, 1))
lab <- c(intercept = "National rate", region = "Regional means", covariates = "Travel time + density",
         main = "Main spatial model", main_region = "Main + region effects", ml_alone_rf = "Random forest alone",
         ml_alone_gbm = "Gradient boosting alone", ml_alone_ridge = "Ridge logistic alone",
         ml_alone_stack = "ML stack alone", ml_glm = "ML stack, no spatial field",
         stacked = "ML stack + spatial field", stacked_cov = "ML stack + covariates + field")
s2 <- summ[order(-summ$mae), ]
col <- ifelse(grepl("ml|stacked", s2$model), "#8e44ad", ifelse(s2$model == "main", "#1b6ca8", "grey60"))
barplot(s2$mae, names.arg = lab[s2$model], horiz = TRUE, las = 1, col = col, border = NA, cex.names = 0.8,
        xlab = "Mean absolute error, held-out district slip rate", xlim = c(0, max(s2$mae) * 1.1),
        main = "Leave-one-district-out: spatial models vs machine learning", cex.main = 0.95)
abline(v = read.csv("data/derived/district_cv_noise.csv")$noise_floor_mae, lty = 2)
mtext("dashed: sampling noise floor", side = 1, line = 3.3, adj = 1, cex = 0.7)
dev.off()
cat("done\n")
