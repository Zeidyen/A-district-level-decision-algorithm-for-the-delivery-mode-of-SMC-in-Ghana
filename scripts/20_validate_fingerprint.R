# External check of the step 2 map: does the DHIMS SMC fingerprint (19_dhims_qa.R; lower =
# stronger SMC signal in routine case data) track the predicted per-contact slip q by
# district (06_inla_timeliness.R)? Independent data (facility case reports vs DHS cards)
# and a different delivery channel (campaign SMC vs routine vaccination).
suppressMessages({library(dplyr)})
norm <- function(x) gsub("[^a-z]", "", tolower(x))
fp <- read.csv("data/derived/dhims_smc_fingerprint.csv")
gs <- read.csv("data/derived/district_gap_slip.csv")          # q by district (u5-weighted)
s4 <- read.csv("data/derived/district_step4.csv")
rb <- read.csv("data/derived/district_robustness.csv")
# DHIMS -> geoBoundaries ADM2 names where spelling differs
alias <- c(gushiegu = "gushegu", tatalesangule = "tatalesanguli", sagnarigu = "sagnerigu",
           tamale = "tamalemetropolitan", yendi = "yendimunicipal", savelugu = "savelugu",
           lambussie = "lambussiekarni", bolgatangaeast = "bolgaeast", kasenanankana = "kasenanankanaeast",
           nadowlikaleo = "nadowlikaleo", westmamprusi = "westmamprusimunicipal",
           eastgonja = "eastgonjamunicipal")
fp$key <- norm(fp$District); fp$key <- ifelse(fp$key %in% names(alias), alias[fp$key], fp$key)
gs$key <- norm(gs$district); s4$key <- norm(s4$district); rb$key <- norm(rb$district)
m <- fp |> inner_join(gs[, c("key", "gap_slip")], by = "key") |>
  inner_join(s4[, c("key", "ratio", "class")], by = "key") |> inner_join(rb[, c("key", "robust")], by = "key")
cat("matched", nrow(m), "of", nrow(fp), "DHIMS districts; unmatched:", paste(setdiff(fp$District, m$District), collapse = ", "), "\n")
stopifnot(nrow(m) == nrow(fp))
ct <- function(a, b) { t <- cor.test(a, b, method = "spearman", exact = FALSE); sprintf("rho %.2f (p = %.3g)", t$estimate, t$p.value) }
cat("Fingerprint vs predicted q:            ", ct(m$index, m$gap_slip), "\n")
cat("Fingerprint vs step 4 routine/campaign:", ct(m$index, m$ratio), "\n")
cat("Within region (rank residuals):        ",
    ct(m$index - ave(m$index, m$Region), m$gap_slip - ave(m$gap_slip, m$Region)), "\n")
cat("\nFingerprint by robustness group (median, n):\n")
print(m |> group_by(robust) |> summarise(n = n(), median_index = round(median(index), 3)) |> as.data.frame(), row.names = FALSE)
cat("\nFingerprint by step 4 class (main):\n")
print(m |> group_by(class) |> summarise(n = n(), median_index = round(median(index), 3)) |> as.data.frame(), row.names = FALSE)
write.csv(m, "data/derived/fingerprint_vs_map.csv", row.names = FALSE)

# The cross-sectional index is confounded by seasonality (19: placebo). The clean test is the
# before/after SMC change in the old Northern region (29 districts) against the map.
w <- read.csv("data/derived/dhims_prepost_oldnorthern.csv") |>
  inner_join(m[, c("District", "gap_slip", "ratio", "class", "robust")], by = "District")
st <- function(a, b) { t <- suppressWarnings(cor.test(a, b, method = "spearman")); c(rho = unname(t$estimate), p = t$p.value) }
on <- w$Region %in% c("Northern", "Savannah", "North East")
v <- rbind(
  cross_all55 = st(m$index, m$gap_slip),
  cross_within_region = st(m$index - ave(m$index, m$Region), m$gap_slip - ave(m$gap_slip, m$Region)),
  prepost_vs_q = st(w$change[on], w$gap_slip[on]),
  prepost_within_region = st((w$change - ave(w$change, w$Region))[on], (w$gap_slip - ave(w$gap_slip, w$Region))[on]),
  prepost_incl_ue = st(w$change, w$gap_slip),
  prepost_incl_ue_within = st(w$change - ave(w$change, w$Region), w$gap_slip - ave(w$gap_slip, w$Region)))
v <- data.frame(test = rownames(v), n = c(nrow(m), nrow(m), sum(on), sum(on), nrow(w), nrow(w)), round(v, 3))
write.csv(v, "data/derived/fingerprint_validation.csv", row.names = FALSE)
print(v, row.names = FALSE)
cat("Old Northern: median change by class:", paste(names(tapply(w$change[on], w$class[on], median)), round(tapply(w$change[on], w$class[on], median), 3), collapse = " | "), "\n")
cat("No SMC drop (change >= 1):", paste(w$District[w$change >= 1], collapse = ", "), "\n")

# ---- under-5-only check (main): drop in the under-5 Jun-Oct share after SMC vs the map ----
u <- read.csv("data/derived/dhims_u5_season_prepost.csv")
mk <- m; mk$adm2_key <- mk$key
u <- inner_join(u, mk[, c("adm2_key", "gap_slip", "ratio", "class")], by = "adm2_key")
on <- u$Region %in% c("Northern", "Savannah", "North East")
# A bigger SMC effect = a more negative change. If the map is right, districts with worse routine
# timeliness (higher q) - which are also those where the campaign struggles - should show less drop:
# expect a POSITIVE correlation between change and q.
vu <- rbind(
  u5_oldnorthern = st(u$change[on], u$gap_slip[on]),
  u5_oldnorthern_within = st((u$change - ave(u$change, u$Region))[on], (u$gap_slip - ave(u$gap_slip, u$Region))[on]),
  u5_incl_ue = st(u$change, u$gap_slip),
  u5_incl_ue_within = st(u$change - ave(u$change, u$Region), u$gap_slip - ave(u$gap_slip, u$Region)))
vu <- data.frame(test = rownames(vu), n = c(sum(on), sum(on), nrow(u), nrow(u)), round(vu, 3))
write.csv(rbind(read.csv("data/derived/fingerprint_validation.csv"), vu), "data/derived/fingerprint_validation.csv", row.names = FALSE)
cat("\nUnder-5-only check (change in Jun-Oct share vs predicted q):\n"); print(vu, row.names = FALSE)
cat("Median change by class (old Northern):", paste(names(tapply(u$change[on], u$class[on], median)),
    sprintf("%+.3f", tapply(u$change[on], u$class[on], median)), collapse = " | "), "\n")
