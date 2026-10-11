# First-order and total-effect Sobol indices of each region's routine-to-campaign ratio.
#
# 21_breakeven_gsa.R gives first-order shares, which sum to less than 1 where inputs interact. Here
# the Jansen (1999) / Saltelli et al. (2010) estimators give, for each input, the first-order
# index S_i (its effect alone) and the total-effect index ST_i (alone plus every interaction);
# ST_i - S_i is the interaction share. Same inputs and distributions as the global sensitivity
# analysis in 21: k (lognormal, median 2, 95% 1.4-2.8), campaign coverage (beta; narrow in Upper
# East), age reach rho ~ U(0.73, 1), routine adherence a_r ~ U(0.744, 0.845), partial-course
# protection f ~ U(0, 1), first-contact reach ~ U(0.9, 1), and the step 2 surface (one of the
# 1,000 posterior draws, drawn uniformly). N = 50,000 base samples (9 x N model evaluations; outputs
# centred before estimation, which leaves the indices unchanged but reduces their variance);
# 95% intervals from 300 bootstrap resamples of the N rows.
#
# Input: data/derived/gsa_region_sums.rds (21). Outputs: data/derived/sobol_total.csv,
# out/sobol_total.png (all five regions)

G <- readRDS("data/derived/gsa_region_sums.rds")
KG <- G$KG; nK <- length(KG); NORTH <- G$NORTH; ND <- dim(G$RS)[2]
K_MEDIAN <- 2; K_SDLOG <- (log(3) - log(1.5)) / 3.92
ADH_C <- 0.845; ADH_R <- ADH_C * 0.8674 / 0.9850
AGE <- read.csv("data/derived/age_reach.csv")
CC_UE <- c(0.87, 0.867, 0.895); CC_OTHER <- c(0.87, 0.75, 0.95)
rbeta_ci <- function(n, m, lo, hi) { v <- ((hi - lo) / 3.92)^2; a <- m * (1 - m) / v - 1; rbeta(n, m * a, (1 - m) * a) }
madh_f <- function(ar, f) (ar + (1 - ar) * f) / (ADH_C + (1 - ADH_C) * f)
IN <- c("k", "cc", "ra", "a_r", "f", "reach", "draw")
LAB <- c(k = "k (routine vs vaccination slip)", cc = "campaign coverage", ra = "age reach (24-59 m)",
         a_r = "routine adherence", f = "partial-course protection", reach = "first-contact reach",
         draw = "step-2 timeliness surface")

sample_inputs <- function(n, r) data.frame(
  k = exp(rnorm(n, log(K_MEDIAN), K_SDLOG)),
  cc = if (r == "Upper East") rbeta_ci(n, CC_UE[1], CC_UE[2], CC_UE[3]) else rbeta_ci(n, CC_OTHER[1], CC_OTHER[2], CC_OTHER[3]),
  ra = runif(n, AGE$lo, 1), a_r = runif(n, ADH_R, ADH_C), f = runif(n, 0, 1), reach = runif(n, 0.9, 1),
  draw = sample.int(ND, n, replace = TRUE))
model <- function(X, r) {
  ri <- match(r, NORTH)
  ki <- pmin(findInterval(X$k, KG), nK - 1); wk <- (X$k - KG[ki]) / (KG[ki + 1] - KG[ki])
  s1 <- G$RS[cbind(ki, X$draw, ri)] * (1 - wk) + G$RS[cbind(ki + 1, X$draw, ri)] * wk
  s2 <- G$RSd[cbind(ki, X$draw, ri)] * (1 - wk) + G$RSd[cbind(ki + 1, X$draw, ri)] * wk
  X$reach * G$P1[r] * madh_f(X$a_r, X$f) / X$cc * (s1 - (1 - X$ra) * s2) / G$RA[r]
}

N <- 50000; B <- 300
set.seed(20261012)
res <- do.call(rbind, lapply(NORTH, function(r) {
  A <- sample_inputs(N, r); Bm <- sample_inputs(N, r)
  fA <- model(A, r); fB <- model(Bm, r)
  fAB <- sapply(IN, function(i) { X <- A; X[[i]] <- Bm[[i]]; model(X, r) })
  m0 <- mean(c(fA, fB)); fA <- fA - m0; fB <- fB - m0; fAB <- fAB - m0   # centring stabilises S
  est <- function(idx) {
    V <- var(c(fA[idx], fB[idx]))
    rbind(S = colMeans(fB[idx] * (fAB[idx, , drop = FALSE] - fA[idx])) / V,          # Saltelli 2010
          ST = colMeans((fA[idx] - fAB[idx, , drop = FALSE])^2) / (2 * V))           # Jansen 1999
  }
  e <- est(seq_len(N))
  bs <- replicate(B, est(sample.int(N, N, replace = TRUE)), simplify = "array")
  data.frame(region = r, input = LAB[IN], S = e["S", ], S_lo = apply(bs["S", , ], 1, quantile, 0.025),
             S_hi = apply(bs["S", , ], 1, quantile, 0.975), ST = e["ST", ],
             ST_lo = apply(bs["ST", , ], 1, quantile, 0.025), ST_hi = apply(bs["ST", , ], 1, quantile, 0.975))
}))
res[, 3:8] <- round(res[, 3:8], 3)
write.csv(res, "data/derived/sobol_total.csv", row.names = FALSE)
cat("First-order (S) and total-effect (ST) indices:\n")
print(res[res$region %in% c("Northern", "Upper East"), c("region", "input", "S", "ST", "ST_lo", "ST_hi")], row.names = FALSE)
cat("\nInteraction share (1 - sum S) by region:",
    paste(NORTH, round(1 - tapply(res$S, res$region, sum)[NORTH], 3), collapse = " | "), "\n")

SHORT <- c("k (routine vs vaccination slip)" = "k", "campaign coverage" = "Campaign coverage",
           "age reach (24-59 m)" = "Age reach (2-4 y)", "routine adherence" = "Routine adherence",
           "partial-course protection" = "Partial-course protection", "first-contact reach" = "First-contact reach",
           "step-2 timeliness surface" = "Timeliness map")
png("out/sobol_total.png", width = 2400, height = 1650, res = 150)
par(mfrow = c(2, 3), mar = c(5, 15, 3.5, 1), cex.axis = 1.35, cex.lab = 1.45)
for (r in NORTH) {
  x <- res[res$region == r, ]; x <- x[order(x$ST), ]
  m <- t(as.matrix(x[, c("S", "ST")])); colnames(m) <- SHORT[x$input]
  bp <- barplot(m, beside = TRUE, horiz = TRUE, las = 1, xlim = c(0, 0.6), col = c("#1b6ca8", "#f0b44c"), border = NA,
                cex.names = 1.45, xlab = "Sobol index", main = r, cex.main = 1.8)
  segments(x$ST_lo, bp[2, ], x$ST_hi, bp[2, ], lwd = 1.6)
}
plot.new()
legend("center", c("First-order index (S)", "Total-effect index (ST)\nwith 95% bootstrap interval"),
       fill = c("#1b6ca8", "#f0b44c"), border = NA, bty = "n", cex = 1.8, y.intersp = 1.6)
dev.off()
cat("done\n")
