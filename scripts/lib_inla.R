# Shared pieces of the step 2 timeliness model, so the main fit (06), the
# checks and the cross-validation (08) use the same data, mesh and model.
# Same model as 06_inla_timeliness.R; see that file for the rationale.

suppressMessages({library(INLA); library(sf); library(terra); library(dplyr)})
sf_use_s2(FALSE)

NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
UTM <- "EPSG:32630"
km_xy <- function(lon, lat)
  st_coordinates(st_transform(st_as_sf(data.frame(lon, lat), coords = 1:2, crs = 4326), UTM)) / 1000

# Outcome: per-contact slip counts (04_timeliness.R): y slipped of n monthly contacts.
# `outcome` picks the threshold column: slip56 (main), slip49, slip42.
load_clusters <- function(outcome = "slip56") {
  cl <- read.csv("data/derived/dhs_clusters.csv") |>
    inner_join(read.csv("data/derived/dhs_clusters_timeliness.csv"), by = "cluster") |>
    filter(contacts > 0) |>
    mutate(y = .data[[outcome]], n = contacts) |>
    mutate(walk = log2(pmax(tt_walk_any, 5)), motor = log2(pmax(tt_motor_any, 5)),
           region = factor(region))
  xy <- km_xy(cl$LONGNUM, cl$LATNUM)
  cl$X <- xy[, 1]; cl$Y <- xy[, 2]
  cl
}

make_mesh <- function(cl) {
  gha <- st_transform(st_union(st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)), UTM)
  bnd <- inla.sp2segment(as(st_sf(geometry = gha / 1000), "Spatial"))
  inla.mesh.2d(loc = cbind(cl$X, cl$Y), boundary = bnd, max.edge = c(25, 100),
               cutoff = 5, offset = c(10, 50))
}

region_dummies <- function(region) {
  lv <- setdiff(sort(unique(as.character(region))), "Greater Accra")
  d <- sapply(lv, function(l) as.integer(region == l))
  colnames(d) <- paste0("rg_", make.names(lv))
  as.data.frame(d)
}

# spec: list(tt = "walk" | "motor", region = TRUE/FALSE, spatial = TRUE/FALSE,
#            covars = TRUE/FALSE, ml = TRUE/FALSE). test: indices whose outcome is hidden.
# prior_range (km, default 30) sets P(range < prior_range) = 0.05 for the Matern field.
# ml = TRUE adds the column cl$ml (a stacked machine-learning prediction on the logit scale,
# 25_ml_stack.py) as a covariate; draw_p() then needs loc$ml too.
fit_model <- function(cl, mesh, spec, test = integer(0)) {
  tr <- setdiff(seq_len(nrow(cl)), test)
  ctr <- c(tt = mean(cl[[spec$tt]][tr]), lpd = mean(cl$logpopdens[tr]),
           ml = if (isTRUE(spec$ml)) mean(cl$ml[tr]) else 0)
  y <- cl$y; y[test] <- NA
  eff <- data.frame(b0 = 1, x1 = cl[[spec$tt]] - ctr["tt"], x2 = cl$logpopdens - ctr["lpd"],
                    clus = seq_len(nrow(cl)))
  # Region effects as dummies against Greater Accra (the intercept), named rg_*.
  if (isTRUE(spec$ml)) eff$x3 <- cl$ml - ctr["ml"]
  rd <- region_dummies(cl$region)
  if (isTRUE(spec$region)) eff <- cbind(eff, rd)
  hyp <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
  terms <- c("-1 + b0",
             if (isTRUE(spec$covars %||% TRUE)) "x1 + x2",
             if (isTRUE(spec$ml)) "x3",
             if (isTRUE(spec$region)) paste(names(rd), collapse = " + "),
             "f(clus, model = 'iid', hyper = hyp)",
             if (isTRUE(spec$spatial)) "f(s, model = spde)")
  f <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
  spde <- inla.spde2.pcmatern(mesh, prior.range = c(spec$prior_range %||% 30, 0.05), prior.sigma = c(1, 0.05))
  if (isTRUE(spec$spatial)) {
    A <- inla.spde.make.A(mesh, loc = cbind(cl$X, cl$Y))
    stk <- inla.stack(data = list(y = y, n = cl$n), A = list(A, 1),
                      effects = list(inla.spde.make.index("s", spde$n.spde), eff), tag = "est")
    dat <- inla.stack.data(stk, spde = spde, hyp = hyp); AA <- inla.stack.A(stk)
    ie <- inla.stack.index(stk, "est")$data
  } else {
    dat <- c(as.list(eff), list(y = y, n = cl$n, hyp = hyp)); AA <- NULL
    ie <- seq_len(nrow(cl))
  }
  fit <- inla(f, family = "binomial", Ntrials = n, data = dat,
              control.predictor = list(A = AA, compute = TRUE, link = 1),
              control.compute = list(config = TRUE, waic = TRUE),
              num.threads = "1:1")   # single thread: reproducible fits
  list(fit = fit, spec = spec, ctr = ctr, mesh = mesh, p = fit$summary.fitted.values$mean[ie])
}

# Prediction grid over the north: covariates, region, under-5 weight.
north_grid <- function() {
  adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
  adm1$region <- sub(" Region$", "", adm1$shapeName)
  nr <- adm1[adm1$region %in% NORTH, ]
  g <- c(rast("data/derived/tt_walk_any_gha.tif"), rast("data/derived/tt_motor_any_gha.tif"),
         rast("data/derived/logpopdens_gha.tif"))
  g <- mask(crop(g, vect(nr)), vect(nr)); names(g) <- c("tw", "tm", "lpd")
  gd <- as.data.frame(g, xy = TRUE, cells = TRUE, na.rm = TRUE)
  gd$region <- nr$region[sapply(st_intersects(st_as_sf(gd, coords = c("x", "y"), crs = 4326), nr),
                                function(i) if (length(i)) i[1] else NA)]
  u5 <- rast("data/raw/gha_t_00_2025_CN_1km_R2025A_UA_v1.tif") +
        rast("data/raw/gha_t_01_2025_CN_1km_R2025A_UA_v1.tif")
  gd$w <- extract(resample(u5, g, method = "bilinear"), gd[, c("x", "y")])[, 2]
  gd$w[is.na(gd$w) | gd$w < 0] <- 0
  gd <- gd[!is.na(gd$region), ]
  gd$walk <- log2(pmax(gd$tw, 5)); gd$motor <- log2(pmax(gd$tm, 5))
  xy <- km_xy(gd$x, gd$y); gd$X <- xy[, 1]; gd$Y <- xy[, 2]
  gd
}

# Posterior draws of p at new locations (each with a fresh cluster effect).
# loc: data frame with X, Y, walk/motor, logpopdens (or lpd), region.
# The iid draws use R's RNG, so callers set.seed(); INLA's sampler takes `seed`.
draw_p <- function(m, loc, n = 300, seed = 20261001L) {
  s <- inla.posterior.sample(n, m$fit, seed = seed)               # INLA ignores set.seed()
  cn <- rownames(s[[1]]$latent)
  x1 <- loc[[m$spec$tt]] - m$ctr["tt"]
  x2 <- (if (!is.null(loc$logpopdens)) loc$logpopdens else loc$lpd) - m$ctr["lpd"]
  A <- if (isTRUE(m$spec$spatial)) inla.spde.make.A(m$mesh, loc = cbind(loc$X, loc$Y))
  is_ <- grep("^s:", cn)
  sapply(s, function(k) {
    lt <- k$latent[, 1]
    eta <- lt["b0:1"] + rnorm(nrow(loc), 0, 1 / sqrt(k$hyperpar[["Precision for clus"]]))
    if (isTRUE(m$spec$covars %||% TRUE)) eta <- eta + lt["x1:1"] * x1 + lt["x2:1"] * x2
    if (isTRUE(m$spec$ml)) eta <- eta + lt["x3:1"] * (loc$ml - m$ctr["ml"])
    if (isTRUE(m$spec$region)) {
      rn <- paste0("rg_", make.names(loc$region), ":1")
      eta <- eta + ifelse(rn %in% names(lt), lt[rn], 0)
    }
    if (isTRUE(m$spec$spatial)) eta <- eta + as.vector(A %*% lt[is_])
    plogis(eta)
  })
}

region_summary <- function(P, region, w) {
  sapply(NORTH, function(r) {
    k <- region == r
    d <- colSums(P[k, , drop = FALSE] * w[k]) / sum(w[k])
    sprintf("%.0f (%.0f-%.0f)", 100 * mean(d), 100 * quantile(d, .025), 100 * quantile(d, .975))
  })
}
