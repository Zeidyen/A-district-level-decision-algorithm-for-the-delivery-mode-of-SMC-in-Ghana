# Step 2: geostatistical model of routine-contact timeliness, predicted across
# the five northern SMC regions.
#
# Outcome: per-contact slip (04_timeliness.R): of the monthly Penta contacts after
# Penta1 (doses 2 and 3), how many slipped by a month or more (gap > 56 days) or
# never happened. Cluster counts y of n contacts; the map is q, the per-contact
# probability, used directly in step 4.
# All Ghana clusters are used for fitting; the map is predicted for the north.
#
#   logit p(s) = b0 + b1 log2(motorised travel time) + b2 log(pop density) + u(s) + e_cluster
#
# Motorised rather than walking time: better WAIC (1307 vs 1314) and spatial CV,
# same regional map (08_checks.R). Walking time is the sensitivity analysis there.
#
# u(s) is a Matern field (SPDE), e an iid cluster effect. Covariates at clusters
# are under-5-weighted means over the DHS displacement buffer (03_dhs_clusters.R). Only
# covariates that exist as rasters can be used, so wealth is not in the model.
#
# Outputs: out/timeliness_map.png, data/derived/gap_slip_{mean,sd,p30}.tif,
# data/derived/district_gap_slip.csv (under-5-weighted, with 95% intervals).

suppressMessages({library(INLA); library(sf); library(terra); library(dplyr)})
sf_use_s2(FALSE)
set.seed(20261001)

NORTH <- c("Northern Region", "Savannah Region", "North East Region",
           "Upper East Region", "Upper West Region")
N_SAMPLES <- 500
UTM <- "EPSG:32630"                      # UTM 30N, kilometres below

# ---- data ---------------------------------------------------------------------
cl <- read.csv("data/derived/dhs_clusters.csv") |>
  inner_join(read.csv("data/derived/dhs_clusters_timeliness.csv"), by = "cluster") |>
  filter(contacts > 0) |>
  mutate(log2_tt = log2(pmax(tt_motor_any, 5)), y = slip56, n = contacts)
xy <- st_coordinates(st_transform(st_as_sf(cl, coords = c("LONGNUM", "LATNUM"), crs = 4326), UTM)) / 1000
cat(nrow(cl), "clusters,", sum(cl$n), "monthly contacts,", sum(cl$y), "slipped\n")

# Centre covariates so the intercept is interpretable and the sampler is stable.
ctr <- c(log2_tt = mean(cl$log2_tt), logpopdens = mean(cl$logpopdens))
cl$x1 <- cl$log2_tt - ctr["log2_tt"]
cl$x2 <- cl$logpopdens - ctr["logpopdens"]

# ---- mesh and SPDE ------------------------------------------------------------
gha <- st_transform(st_union(st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)), UTM)
bnd <- inla.sp2segment(as(st_sf(geometry = gha / 1000), "Spatial"))
mesh <- inla.mesh.2d(loc = xy, boundary = bnd, max.edge = c(25, 100), cutoff = 5, offset = c(10, 50))
cat("mesh nodes:", mesh$n, "\n")
# PC priors: P(range < 30 km) = 0.05; P(sd > 1) = 0.05.
spde <- inla.spde2.pcmatern(mesh, prior.range = c(30, 0.05), prior.sigma = c(1, 0.05))
A <- inla.spde.make.A(mesh, loc = xy)
idx <- inla.spde.make.index("s", spde$n.spde)

stk <- inla.stack(data = list(y = cl$y, n = cl$n), A = list(A, 1),
                  effects = list(idx, data.frame(b0 = 1, x1 = cl$x1, x2 = cl$x2,
                                                 clus = seq_len(nrow(cl)))),
                  tag = "est")
f <- y ~ -1 + b0 + x1 + x2 + f(s, model = spde) +
  f(clus, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05))))
fit <- inla(f, family = "binomial", Ntrials = n, data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
            control.compute = list(config = TRUE, waic = TRUE, cpo = TRUE),
            num.threads = "1:1")   # single thread: reproducible fits

fx <- fit$summary.fixed
cat("\nFixed effects (odds ratios, 95% CrI):\n")
cat(sprintf("  per doubling of travel time: %.2f (%.2f-%.2f)\n",
            exp(fx["x1", "mean"]), exp(fx["x1", "0.025quant"]), exp(fx["x1", "0.975quant"])))
cat(sprintf("  per unit log pop density:     %.2f (%.2f-%.2f)\n",
            exp(fx["x2", "mean"]), exp(fx["x2", "0.025quant"]), exp(fx["x2", "0.975quant"])))
hp <- fit$summary.hyperpar
print(round(hp[, c("mean", "0.025quant", "0.975quant")], 2))
cat("WAIC:", round(fit$waic$waic, 1), "| CPO failures:", sum(fit$cpo$failure > 0, na.rm = TRUE), "\n")

# Covariate-only comparison: does the spatial field earn its place?
fit0 <- inla(y ~ 1 + x1 + x2 + f(clus, model = "iid",
               hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05)))),
             family = "binomial", Ntrials = n,
             data = transform(cl, clus = seq_len(nrow(cl))),
             control.compute = list(waic = TRUE),
             num.threads = "1:1")   # single thread: reproducible fits
cat("WAIC without spatial field:", round(fit0$waic$waic, 1), "\n")
write.csv(data.frame(term = c("travel_time_doubling", "log_pop_density"),
                     or = exp(fx[c("x1", "x2"), "mean"]), lo = exp(fx[c("x1", "x2"), "0.025quant"]),
                     hi = exp(fx[c("x1", "x2"), "0.975quant"]),
                     waic = fit$waic$waic, waic_no_field = fit0$waic$waic),
          "data/derived/step2_fixed.csv", row.names = FALSE)

# In-sample check: fitted vs observed at the clusters, child-weighted, by region.
ie <- inla.stack.index(stk, "est")$data
cl$fitted <- fit$summary.fitted.values$mean[ie]
cat("\nObserved vs fitted at clusters (child-weighted %):\n")
print(cl |> group_by(region) |>
        summarise(observed = round(100 * sum(y) / sum(n)),
                  fitted = round(100 * sum(fitted * n) / sum(n))) |>
        filter(region %in% sub(" Region", "", NORTH)) |> as.data.frame(), row.names = FALSE)

# ---- prediction over the north ------------------------------------------------
north <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE) |> filter(shapeName %in% NORTH)
tt <- rast("data/derived/tt_motor_any_gha.tif")
lpd <- rast("data/derived/logpopdens_gha.tif")
grid <- mask(crop(c(tt, lpd), vect(north)), vect(north))
names(grid) <- c("tt", "lpd")
gd <- as.data.frame(grid, xy = TRUE, cells = TRUE, na.rm = TRUE)
gxy <- st_coordinates(st_transform(st_as_sf(gd, coords = c("x", "y"), crs = 4326), UTM)) / 1000
Ap <- inla.spde.make.A(mesh, loc = gxy)
X1 <- log2(pmax(gd$tt, 5)) - ctr["log2_tt"]
X2 <- gd$lpd - ctr["logpopdens"]
cat("prediction pixels:", nrow(gd), "\n")

# Joint posterior samples. Each draw adds a fresh cluster-level iid effect, so
# the map describes a new community, not the average of communities.
samp <- inla.posterior.sample(N_SAMPLES, fit, seed = 20261001L)  # INLA ignores set.seed()
cn <- rownames(samp[[1]]$latent)
i_s <- grep("^s:", cn); i_b0 <- match("b0:1", cn); i_x1 <- match("x1:1", cn); i_x2 <- match("x2:1", cn)
P <- matrix(NA_real_, nrow(gd), N_SAMPLES)
for (k in seq_len(N_SAMPLES)) {
  lt <- samp[[k]]$latent
  tau <- samp[[k]]$hyperpar[["Precision for clus"]]
  eta <- lt[i_b0] + lt[i_x1] * X1 + lt[i_x2] * X2 + as.vector(Ap %*% lt[i_s]) +
         rnorm(nrow(gd), 0, 1 / sqrt(tau))
  P[, k] <- plogis(eta)
}

out <- rast(grid[[1]]); names(out) <- "x"
put <- function(v, nm) { r <- out; r[gd$cell] <- v; names(r) <- nm; r }
r_mean <- put(rowMeans(P), "gap_slip_mean")
r_sd   <- put(apply(P, 1, sd), "gap_slip_sd")
r_p30  <- put(rowMeans(P > 0.30), "p_gap_slip_gt30")
for (r in list(r_mean, r_sd, r_p30))
  writeRaster(r, sprintf("data/derived/%s.tif", names(r)), overwrite = TRUE)

# ---- district summaries, weighted by children under 5 -------------------------
u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") +
      rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
w <- extract(resample(u5, grid, method = "bilinear"), gd[, c("x", "y")])[, 2]
w[is.na(w) | w < 0] <- 0

adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)]
adm2 <- filter(adm2, region %in% NORTH)
dist_id <- sapply(st_intersects(st_as_sf(gd, coords = c("x", "y"), crs = 4326), adm2),
                  function(i) if (length(i)) i[1] else NA)

ds <- lapply(seq_len(nrow(adm2)), function(j) {
  k <- which(dist_id == j)
  draws <- colSums(P[k, , drop = FALSE] * w[k]) / sum(w[k])
  data.frame(region = adm2$region[j], district = adm2$shapeName[j], u5 = round(sum(w[k])),
             gap_slip = round(mean(draws), 3),
             lo = round(quantile(draws, 0.025), 3), hi = round(quantile(draws, 0.975), 3))
})
ds <- do.call(rbind, ds) |> arrange(desc(gap_slip))
write.csv(ds, "data/derived/district_gap_slip.csv", row.names = FALSE)

reg <- lapply(NORTH, function(rg) {
  k <- which(adm2$region[dist_id] == rg)
  draws <- colSums(P[k, ] * w[k]) / sum(w[k])
  sprintf("%-18s %.0f%% (%.0f-%.0f)", rg, 100 * mean(draws),
          100 * quantile(draws, 0.025), 100 * quantile(draws, 0.975))
})
cat("\nPredicted per-contact slip probability q (u5-weighted):\n ",
    paste(unlist(reg), collapse = "\n  "), "\n")
cat("\nHighest districts:\n"); print(head(ds, 8), row.names = FALSE)
cat("Lowest districts:\n");  print(tail(ds, 5), row.names = FALSE)

# ---- map ----------------------------------------------------------------------
png("out/timeliness_map.png", width = 2400, height = 1000, res = 150)
par(mfrow = c(1, 3), mar = c(1, 1, 3, 6))
nl <- st_geometry(north); d2 <- st_geometry(adm2)
cpts <- st_as_sf(cl, coords = c("LONGNUM", "LATNUM"), crs = 4326)
cpts <- cpts[st_intersects(cpts, st_union(north), sparse = FALSE)[, 1], ]
plot(clamp(r_mean, 0, 0.7), col = hcl.colors(50, "YlOrRd", rev = TRUE), range = c(0, 0.7),  # values > 0.7 shown at the top colour
     main = "Predicted chance a monthly contact slips (q)", axes = FALSE, box = FALSE)
plot(d2, add = TRUE, border = "grey50", lwd = 0.3); plot(nl, add = TRUE, lwd = 1)
plot(st_geometry(cpts), add = TRUE, pch = 3, cex = 0.4)
plot(r_sd, col = hcl.colors(50, "Purples", rev = TRUE),
     main = "Posterior SD", axes = FALSE, box = FALSE)
plot(d2, add = TRUE, border = "grey50", lwd = 0.3); plot(nl, add = TRUE, lwd = 1)
plot(r_p30, col = hcl.colors(50, "RdYlBu", rev = TRUE), range = c(0, 1),
     main = "P(q > 30%)", axes = FALSE, box = FALSE)
plot(d2, add = TRUE, border = "grey50", lwd = 0.3); plot(nl, add = TRUE, lwd = 1)
dev.off()
cat("done\n")
