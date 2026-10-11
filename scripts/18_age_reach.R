# Does routine contact fall off for older under-5s? SMC covers 3-59 months, but routine
# immunisation contacts end at ~18 months (MR2). Vitamin A in the last 6 months (h34) is
# given at child welfare clinics (and in child health campaigns), and deworming (h43) at
# CWC and schools: the nearest DHS measures of routine contact for children past infancy.
# Northern regions, living children, survey-weighted.
# Output: data/derived/age_reach.csv (ratio older / younger, used by 10_step4_impact.R "age").
suppressMessages({library(haven); library(dplyr)})
north <- c("northern", "savannah", "north east", "upper east", "upper west")
k <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA") |>
  mutate(age_m = v008 - b3, region = as.character(as_factor(v024)), w = v005 / 1e6) |>
  filter(b5 == 1, region %in% north, age_m >= 6, age_m <= 59) |>
  mutate(band = cut(age_m, c(6, 12, 24, 36, 48, 60), right = FALSE,
                    labels = c("6-11", "12-23", "24-35", "36-47", "48-59")),
         vita = as.integer(h34 == 1), deworm = as.integer(h43 == 1))
for (x in c("h34", "h43")) cat(x, ":", attr(k[[x]], "label"), "| codes", paste(sort(unique(k[[x]])), collapse = ","), "\n")
wm <- function(y, w) { i <- !is.na(y); sum(y[i] * w[i]) / sum(w[i]) }
tab <- k |> group_by(band) |> summarise(n = n(), vitA = round(100 * wm(vita, w), 1),
                                        deworm = round(100 * wm(deworm, w), 1), .groups = "drop")
print(as.data.frame(tab), row.names = FALSE)
write.csv(tab, "data/derived/age_bands_vitA.csv", row.names = FALSE)
by_reg <- k |> mutate(old = age_m >= 24) |> group_by(region, old) |>
  summarise(vitA = wm(vita, w), .groups = "drop") |>
  tidyr::pivot_wider(names_from = old, values_from = vitA, names_prefix = "old_")
by_reg$ratio <- round(by_reg$old_TRUE / by_reg$old_FALSE, 3)
print(as.data.frame(by_reg), row.names = FALSE)

# Ratio used in step 4: contact at 24-59 m relative to 12-23 m (pooled north), with a
# cluster bootstrap interval.
ratio <- function(d) wm(d$vita[d$age_m >= 24], d$w[d$age_m >= 24]) / wm(d$vita[d$age_m %in% 12:23], d$w[d$age_m %in% 12:23])
kk <- k |> select(v001, age_m, vita, w)                 # bootstrap on the 4 columns needed
byc <- split(kk, kk$v001)
set.seed(1)
bs <- replicate(500, ratio(do.call(rbind, byc[sample(names(byc), replace = TRUE)])))
r <- ratio(k); ci <- quantile(bs, c(.025, .975))
cat(sprintf("\nVitamin A, 24-59 m relative to 12-23 m (north): %.3f (95%% CI %.3f-%.3f)\n", r, ci[1], ci[2]))
rd <- wm(k$deworm[k$age_m >= 24], k$w[k$age_m >= 24]) / wm(k$deworm[k$age_m %in% 12:23], k$w[k$age_m %in% 12:23])
cat(sprintf("Deworming, same ratio: %.3f\n", rd))
write.csv(data.frame(measure = "vitA_24_59_vs_12_23", ratio = round(r, 3), lo = round(ci[1], 3), hi = round(ci[2], 3)),
          "data/derived/age_reach.csv", row.names = FALSE)
