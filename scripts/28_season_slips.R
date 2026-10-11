# Do routine contacts slip more in the SMC season (June-September)?
#
# The step 2 model uses Penta contacts due at any time of year, but SMC runs only June-September
# (rainy and farming season). If contacts slip more then, routine SMC would do worse than step 4
# assumes. Card dates give the calendar month each contact fell due (27a_slip_features.R).
#
# Model (northern contacts; all Ghana as a check): mixed logistic regression
#   slip ~ due in Jun-Sep + third dose + log2 walking time + wealth + mother's education
#          + region + year due + (1 | cluster)
# The year term absorbs the 2020 COVID disruption, which fell in particular months. A second
# model replaces the season indicator with the calendar month (reference: October-May pooled).
#
# Scale factor for step 4: f = mean predicted slip if every northern contact fell due in Jun-Sep
# / mean predicted slip as observed (marginal standardisation), with 1,000 draws from the fixed
# effects' sampling distribution. f is about 1 (no season effect), so no step 4 scenario uses it.
#
# Outputs: data/derived/season_slip_or.csv, season_slip_factor.csv, season_slip_months.csv,
#          out/season_slips.png

suppressMessages({library(glmmTMB); library(dplyr); library(MASS)})
d <- read.csv("data/derived/slip_contacts_features.csv") |>
  mutate(contact3 = as.integer(contact == 3), log2_tt = log2(pmax(tt_walk_any, 5)),
         region = factor(region), year = factor(due_year))
n <- filter(d, north) |> mutate(region = droplevels(region))
fit <- function(df, f) glmmTMB(f, family = binomial, data = df)
F1 <- slip ~ smc_season + contact3 + log2_tt + wealth + mother_edu + region + year + (1 | cluster)
m_n <- fit(n, F1)
m_all <- fit(d, F1)
or <- function(m, term, lab) {
  s <- summary(m)$coefficients$cond
  data.frame(model = lab, term = term, or = exp(s[term, 1]),
             lo = exp(s[term, 1] - 1.96 * s[term, 2]), hi = exp(s[term, 1] + 1.96 * s[term, 2]),
             p = s[term, 4])
}
ors <- rbind(or(m_n, "smc_season", "north, adjusted"), or(m_all, "smc_season", "all Ghana, adjusted"),
             or(fit(n, slip ~ smc_season + (1 | cluster)), "smc_season", "north, unadjusted"))
ors[, 3:6] <- round(ors[, 3:6], 3)
write.csv(ors, "data/derived/season_slip_or.csv", row.names = FALSE)
cat("Odds ratio of a slip, contact due in June-September vs other months:\n"); print(ors, row.names = FALSE)

# Month profile (north): June..September each against October-May pooled.
n$mon <- factor(ifelse(n$due_month %in% 6:9, month.abb[n$due_month], "Oct-May"), c("Oct-May", month.abb[6:9]))
m_m <- fit(n, slip ~ mon + contact3 + log2_tt + wealth + mother_edu + region + year + (1 | cluster))
sm <- summary(m_m)$coefficients$cond
mo <- data.frame(month = month.abb[6:9], or = exp(sm[paste0("mon", month.abb[6:9]), 1]),
                 lo = exp(sm[paste0("mon", month.abb[6:9]), 1] - 1.96 * sm[paste0("mon", month.abb[6:9]), 2]),
                 hi = exp(sm[paste0("mon", month.abb[6:9]), 1] + 1.96 * sm[paste0("mon", month.abb[6:9]), 2]),
                 contacts = as.vector(table(n$mon)[month.abb[6:9]]))
mo[, 2:4] <- round(mo[, 2:4], 3)
write.csv(mo, "data/derived/season_slip_months.csv", row.names = FALSE)
cat("\nBy month (vs October-May), north:\n"); print(mo, row.names = FALSE)

# Scale factor f for step 4, by marginal standardisation over the fixed effects (cluster effects
# at zero), with sampling uncertainty from the fixed-effect covariance.
set.seed(20261011)
X <- getME(m_n, "X")                                          # fixed-effect design, rows as fitted
X1 <- X; X1[, "smc_season"] <- 1
b <- fixef(m_n)$cond; V <- vcov(m_n)$cond
B <- mvrnorm(1000, b, V)
f <- apply(B, 1, function(bb) mean(plogis(X1 %*% bb)) / mean(plogis(X %*% bb)))
fs <- data.frame(f_median = median(f), f_lo = quantile(f, 0.025), f_hi = quantile(f, 0.975),
                 f_sdlog = sd(log(f)))
fs[] <- round(fs, 4)
write.csv(fs, "data/derived/season_slip_factor.csv", row.names = FALSE)
cat("\nScale factor for q if every contact fell in the SMC season (north):", fs$f_median,
    sprintf("(95%% %.3f-%.3f)\n", fs$f_lo, fs$f_hi))

# Figure: raw slip share by month due (north), with binomial 95% CI and the SMC season shaded.
r <- n |> group_by(due_month) |> summarise(k = sum(slip), m = n(), .groups = "drop") |>
  mutate(p = k / m, lo = qbeta(0.025, k + 0.5, m - k + 0.5), hi = qbeta(0.975, k + 0.5, m - k + 0.5))
png("out/season_slips.png", width = 1700, height = 950, res = 170)
par(mar = c(4, 4.5, 3, 1))
plot(NA, xlim = c(0.5, 12.5), ylim = c(0, max(r$hi) * 1.05), xaxt = "n", las = 1, xlab = "Month the dose fell due",
     ylab = "Share of routine contacts that slipped",
     main = sprintf("Northern routine contacts by month due (adjusted OR, June-September: %.2f, 95%% CI %.2f-%.2f)",
                    ors$or[1], ors$lo[1], ors$hi[1]), cex.main = 0.85)
rect(5.5, 0, 9.5, max(r$hi) * 1.05, col = "grey92", border = NA)
text(7.5, max(r$hi) * 1.0, "SMC season (Jun-Sep)", cex = 0.8, col = "grey30")
axis(1, 1:12, month.abb)
segments(r$due_month, r$lo, r$due_month, r$hi, col = "#1b6ca8", lwd = 2)
points(r$due_month, r$p, pch = 16, col = "#1b6ca8", cex = sqrt(r$m) / 8)
abline(h = mean(n$slip), lty = 2, col = "grey40")
dev.off()
cat("done\n")
