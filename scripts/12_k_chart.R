# Summary chart: district classes under each step 4 variant (data/derived/step4_sensitivity.csv,
# from 14_sensitivity.R). Never typed numbers.
s <- read.csv("data/derived/step4_sensitivity.csv")
s <- s[nrow(s):1, ]                                    # main at the top
m <- t(as.matrix(s[, c("green", "amber", "red")])); colnames(m) <- s$label
cols <- c("#5cb85c", "#f6c35b", "#d9534f")
png("out/summary_k_chart.png", width = 2000, height = 1250, res = 200)
par(mar = c(6.5, 19, 3.5, 1))
bp <- barplot(m, horiz = TRUE, col = cols, border = NA, las = 1, cex.names = 0.72, xlim = c(0, sum(m[, 1])),
              main = sprintf("Green districts: %d in the main analysis, %d-%d across all %d variants",
                             s$green[s$scenario == "main"], min(s$green), max(s$green), nrow(s)), cex.main = 0.95)
cum <- apply(m, 2, cumsum)
for (j in seq_len(ncol(m))) for (i in 1:3) if (m[i, j] >= 2) text(cum[i, j] - m[i, j] / 2, bp[j], m[i, j], cex = 0.7)
mtext(sprintf("Number of the %d northern districts", sum(m[, 1])), side = 1, line = 2.3, cex = 0.8)
par(xpd = TRUE)
legend("bottom", c("Green: routine averts >= 90% of campaign cases", "Amber: 70-90%", "Red: < 70%"),
       fill = cols, border = NA, bty = "n", horiz = TRUE, cex = 0.7, inset = c(0, -0.2))
dev.off()
