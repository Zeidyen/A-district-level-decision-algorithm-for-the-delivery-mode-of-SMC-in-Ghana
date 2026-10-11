# DHIMS2 district-month malaria data (all Ghana, 2015-2025; 00_dhims_prepare.R): quality
# check for the five northern SMC regions, and the SMC fingerprint.
#
# SMC fingerprint: SMC (doses June-September, each protecting ~28 days) protects 3-59 m in Jun-Oct; 5-14 year-olds share the same
# transmission but get no SMC. Where SMC works, the under-5 : 5-14 ratio of confirmed cases
# drops in Jun-Oct relative to the other months. Index = ratio(Jun-Oct) / ratio(other
# months), per district, 2020-24 pooled. Lower = stronger SMC signal.
# The index is NOT a clean SMC measure. A placebo with two age groups that never get SMC
# (5-14 : 15+, same districts and months) is ~1.37: in peak season, 5-14-year-olds' cases
# surge, which by itself lowers the u5 : 5-14 ratio. Before SMC began (old Northern region,
# 2015-18) the index was already ~0.89.
# SMC start years (confirmed by Zeid, 9 Oct 2026): Upper West 2015; Upper East 2016;
# Northern 2019, including Savannah and North East (split from Northern in 2019). So most of the cross-sectional index reflects seasonality; the SMC effect is the
# before/after change in the same districts. Southern districts are NOT a control: they do
# not qualify for SMC because transmission there is perennial (Zeid, 9 Oct 2026).
# Outputs: data/derived/dhims_smc_fingerprint.csv (55 northern districts),
#          data/derived/dhims_fingerprint_by_year.csv (all districts x year),
#          data/derived/dhims_region_stats_all.csv
suppressMessages({library(dplyr); library(tidyr)})
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
CONTROL <- character(0)   # no valid non-SMC control region (see header)
da <- read.csv("data/raw/dhims/dhims_national_2015_2025.csv", check.names = FALSE)
d <- filter(da, Region %in% NORTH)
cat("North:", nrow(d), "district-months,", n_distinct(d$District), "districts;",
    "zero conf_tot:", sum(d$conf_tot == 0), "\n")

cat("\nConfirmed cases (thousands) by region and year:\n")
print(d |> group_by(Region, year) |> summarise(c = round(sum(conf_tot) / 1000), .groups = "drop") |>
        pivot_wider(names_from = year, values_from = c) |> as.data.frame(), row.names = FALSE)

reg <- da |> filter(year %in% 2020:2024) |> group_by(Region) |>
  summarise(u5_share = sum(conf_u5) / sum(conf_tot), u5_pop_share = sum(Pop_5) / sum(`Pop (total)`),
            season_all = sum(conf_tot[m %in% 6:10]) / sum(conf_tot),
            season_u5 = sum(conf_u5[m %in% 6:10]) / sum(conf_u5),
            tpr_u5 = sum(conf_u5) / sum(test_u5), test_rate_u5 = sum(test_u5) / sum(susp_u5),
            u5_inc = sum(conf_u5) / (sum(Pop_5) / 12) / 5, .groups = "drop") |>
  mutate(group = ifelse(Region %in% NORTH, "north (SMC)", ifelse(Region %in% CONTROL, "control (no SMC)", "south (perennial)")))
write.csv(reg, "data/derived/dhims_region_stats_all.csv", row.names = FALSE)
cat("\nBy region, 2020-24:\n"); print(reg |> mutate(across(where(is.numeric), ~ round(.x, 3))) |> as.data.frame(), row.names = FALSE)

# ---- fingerprint ---------------------------------------------------------------------
fpi <- function(x) x |> mutate(smc = m %in% 6:10) |>
  group_by(Region, District, smc, .add = TRUE) |> summarise(u5 = sum(conf_u5), o = sum(conf_5_14), .groups = "drop_last") |>
  summarise(index = (u5[smc] / o[smc]) / (u5[!smc] / o[!smc]), .groups = "drop")
fp_all <- fpi(da |> filter(year %in% 2020:2024) |> group_by(Region, District))
fp_all$group <- ifelse(fp_all$Region %in% NORTH, "north", ifelse(fp_all$Region %in% CONTROL, "control", "south"))
cat("\nFingerprint, 2020-24, by region:\n")
print(fp_all |> filter(group == "north") |> group_by(Region) |>
        summarise(n = n(), median = round(median(index), 3), min = round(min(index), 3), max = round(max(index), 3)) |>
        as.data.frame(), row.names = FALSE)
fp <- fp_all |> filter(group == "north") |>
  mutate(index = round(index, 3)) |> arrange(index)
write.csv(fp |> select(Region, District, index), "data/derived/dhims_smc_fingerprint.csv", row.names = FALSE)
cat("\nNorthern districts (index; lower = stronger SMC signal):\n")
print(as.data.frame(fp[, c("Region", "District", "index")]), row.names = FALSE)

# By year (district and region), for the pre/post analysis once SMC start years are confirmed.
by_year <- da |> group_by(Region, District, year) |> fpi() |>
  mutate(group = ifelse(Region %in% NORTH, "north", ifelse(Region %in% CONTROL, "control", "south")))
write.csv(by_year, "data/derived/dhims_fingerprint_by_year.csv", row.names = FALSE)
cat("\nMedian district index by year, northern regions:\n")
print(by_year |> filter(group == "north") |>
        group_by(Region, year) |> summarise(i = round(median(index), 2), .groups = "drop") |>
        pivot_wider(names_from = year, values_from = i) |> as.data.frame(), row.names = FALSE)

# ---- placebo and before/after ----------------------------------------------------------
ix <- function(a, b, smc) (sum(a[smc]) / sum(b[smc])) / (sum(a[!smc]) / sum(b[!smc]))
per <- function(y) ifelse(y %in% 2015:2018, "pre", ifelse(y %in% 2020:2024, "post", NA))
pp <- d |> mutate(period = per(year), smc = m %in% 6:10) |> filter(!is.na(period)) |>
  group_by(Region, District, period) |>
  summarise(smc_index = ix(conf_u5, conf_5_14, smc), placebo = ix(conf_5_14, conf_ov15, smc), .groups = "drop")
cat("\nPlacebo (5-14 : 15+) and SMC index, median by region and period (pre 2015-18, post 2020-24):\n")
print(pp |> group_by(Region, period) |> summarise(smc_index = round(median(smc_index), 3), placebo = round(median(placebo), 3), .groups = "drop") |>
        as.data.frame(), row.names = FALSE)
# Before/after where the data include pre-SMC years: pre = 2015 to (start - 1), post = 2020-24.
# Old Northern (Northern, Savannah, North East): 4 pre years. Upper East: 1 (2015), a thin
# baseline. Upper West has none (SMC from 2015, the first year of the data).
START <- c(`Upper West` = 2015, `Upper East` = 2016, Northern = 2019, Savannah = 2019, `North East` = 2019)
w <- d |> mutate(start = START[Region], period = ifelse(year < start, "pre", ifelse(year %in% 2020:2024, "post", NA)),
                 smc = m %in% 6:10) |> filter(!is.na(period)) |>
  group_by(Region, District, period) |>
  summarise(smc_index = ix(conf_u5, conf_5_14, smc), placebo = ix(conf_5_14, conf_ov15, smc),
            pre_years = n_distinct(year), .groups = "drop") |>
  select(-pre_years) |>
  pivot_wider(names_from = period, values_from = c(smc_index, placebo)) |>
  filter(!is.na(smc_index_pre)) |> mutate(change = smc_index_post / smc_index_pre)
on <- w$Region %in% c("Northern", "Savannah", "North East")
cat(sprintf("\nOld Northern (%d districts, 4 pre years): SMC index %.3f -> %.3f (fell in %d); placebo %.3f -> %.3f; change median %.3f\n",
            sum(on), median(w$smc_index_pre[on]), median(w$smc_index_post[on]), sum(w$change[on] < 1),
            median(w$placebo_pre[on]), median(w$placebo_post[on]), median(w$change[on])))
cat(sprintf("Upper East (%d districts, 1 pre year): SMC index %.3f -> %.3f (fell in %d); change median %.3f\n",
            sum(!on), median(w$smc_index_pre[!on]), median(w$smc_index_post[!on]), sum(w$change[!on] < 1), median(w$change[!on])))
write.csv(w, "data/derived/dhims_prepost_oldnorthern.csv", row.names = FALSE)
write.csv(pp, "data/derived/dhims_placebo.csv", row.names = FALSE)

# ---- under-5-only SMC check (main) ----------------------------------------------------
# Uses under-5 cases only: SMC protects under-5s in Jun-Oct, so once it starts, the share of
# a year's under-5 confirmed cases that fall in Jun-Oct should drop. Same districts, before
# (2015 to start - 1) vs after (2020-24). Caveat: rainfall seasonality also varies between
# years, and nothing in under-5 data alone separates that from SMC.
u5s <- d |> mutate(start = START[Region], period = ifelse(year < start, "pre", ifelse(year %in% 2020:2024, "post", NA))) |>
  filter(!is.na(period)) |> group_by(Region, District, adm2_key, period) |>
  summarise(share = sum(conf_u5[m %in% 6:10]) / sum(conf_u5), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = share, names_prefix = "share_") |>
  filter(!is.na(share_pre)) |> mutate(change = share_post - share_pre)
onu <- u5s$Region %in% c("Northern", "Savannah", "North East")
cat(sprintf("\nUnder-5 Jun-Oct share of cases. Old Northern (%d districts, 4 pre years): %.3f -> %.3f, change median %+.3f, fell in %d\n",
            sum(onu), median(u5s$share_pre[onu]), median(u5s$share_post[onu]), median(u5s$change[onu]), sum(u5s$change[onu] < 0)))
cat(sprintf("Upper East (%d districts, 1 pre year): %.3f -> %.3f, change median %+.3f, fell in %d\n",
            sum(!onu), median(u5s$share_pre[!onu]), median(u5s$share_post[!onu]), median(u5s$change[!onu]), sum(u5s$change[!onu] < 0)))
write.csv(u5s, "data/derived/dhims_u5_season_prepost.csv", row.names = FALSE)
