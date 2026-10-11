# Hybrid delivery: how many children could move to routine delivery for how little protection?
#
# Districts are moved from the campaign to routine delivery one at a time, best routine-to-campaign
# ratio first (main analysis, 10_step4_impact.R). After each move: share of eligible children
# moved off the campaign, and cases averted lost as a share of what the all-campaign strategy
# averts. The same order is then re-scored with the age-gap + adherence variant, to show what the
# chosen districts would cost if those assumptions hold.
#
# Reported: the move sets defined by the green line (ratio >= 0.9), by P(green) >= 0.5, and the
# most children that can be moved while losing at most 1%, 2% and 5% of cases averted.
#
# Outputs: data/derived/hybrid_curve.csv, hybrid_summary.csv, out/hybrid_curve.png

suppressMessages(library(dplyr))
m <- read.csv("data/derived/district_step4.csv")
v <- read.csv("data/derived/district_step4_age_adh.csv")
m <- m[order(-m$ratio), ]
m$loss_real <- v$shortfall[match(m$district, v$district)]
tot <- sum(m$averted_campaign); tot_real <- sum(v$averted_campaign); kids <- sum(m$children_3_59m)
curve <- data.frame(step = 0:nrow(m), district = c("", trimws(m$district)), region = c("", m$region),
                    ratio = c(NA, m$ratio),
                    kids_moved = c(0, cumsum(m$children_3_59m)) / kids,
                    loss = c(0, cumsum(m$shortfall)) / tot,
                    loss_age_adh = c(0, cumsum(m$loss_real)) / tot_real)
write.csv(transform(curve, kids_moved = round(kids_moved, 4), loss = round(loss, 4),
                    loss_age_adh = round(loss_age_adh, 4)), "data/derived/hybrid_curve.csv", row.names = FALSE)

pick <- function(lab, n) {
  s <- m[seq_len(n), ]
  data.frame(set = lab, districts = n, children = sum(s$children_3_59m), pct_children = 100 * sum(s$children_3_59m) / kids,
             pct_loss = 100 * sum(s$shortfall) / tot, pct_loss_age_adh = 100 * sum(s$loss_real) / tot_real,
             regions = paste(sort(unique(s$region)), collapse = ", "))
}
budget <- function(b) max(which(curve$loss <= b)) - 1
summ <- rbind(pick("green (ratio >= 0.9)", sum(m$ratio >= 0.9)),
              pick("P(green) >= 0.5", sum(cumsum(m$p_green >= 0.5) == seq_len(nrow(m)) & m$p_green >= 0.5)),
              pick("loss <= 1% of cases averted", budget(0.01)),
              pick("loss <= 2%", budget(0.02)),
              pick("loss <= 5%", budget(0.05)),
              pick("all districts", nrow(m)))
summ[, c("pct_children", "pct_loss", "pct_loss_age_adh")] <- round(summ[, c("pct_children", "pct_loss", "pct_loss_age_adh")], 1)
write.csv(summ, "data/derived/hybrid_summary.csv", row.names = FALSE)
print(summ, row.names = FALSE)
cat("\nFirst districts moved:", paste(head(trimws(m$district), 8), collapse = ", "), "\n")

png("out/hybrid_curve.png", width = 1700, height = 1000, res = 170)
par(mar = c(4.5, 4.5, 3, 1))
plot(100 * curve$kids_moved, 100 * curve$loss, type = "s", lwd = 2.2, col = "#1b6ca8", las = 1,
     xlab = "Eligible children moved from the campaign to routine delivery (%)",
     ylab = "Cases averted lost (% of the all-campaign total)",
     ylim = c(0, 100 * max(curve$loss_age_adh)), main = "Hybrid delivery: districts moved best-first", cex.main = 0.95)
lines(100 * curve$kids_moved, 100 * curve$loss_age_adh, type = "s", lwd = 2.2, col = "#c0392b", lty = 2)
for (b in c(0.01, 0.02, 0.05)) abline(h = 100 * b, col = "grey80", lty = 3)
ng <- sum(m$ratio >= 0.9)
points(100 * curve$kids_moved[ng + 1], 100 * curve$loss[ng + 1], pch = 19, col = "#2e8b57", cex = 1.3)
text(100 * curve$kids_moved[ng + 1], 100 * curve$loss[ng + 1], "all green districts moved", pos = 4, cex = 0.75, col = "#2e8b57")
legend("topleft", c("main analysis", "same districts, with the age gap and adherence loss"),
       col = c("#1b6ca8", "#c0392b"), lty = c(1, 2), lwd = 2.2, bty = "n", cex = 0.8)
dev.off()
cat("done\n")
