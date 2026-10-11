# Contact-level dataset for the "why do routine contacts slip?" analysis (27_slip_drivers.py).
#
# Same contacts and outcome as 04_timeliness.R: children 12-35 months with a card date for
# Penta1; Penta doses 2 and 3 are each a monthly contact, which slips if the gap from the previous
# dose exceeds 56 days or the dose was never given. One row per contact, all Ghana.
#
# Missing-value codes are set per variable (8 and 9 are valid for birth order, household size
# and antenatal visits).
# Features: child (sex, birth order, month of birth, contact number), mother (age, education,
# work, health insurance, media exposure, problems getting care), pregnancy and birth care
# (antenatal visits, place of delivery), household (wealth, size, head's sex, under-5s), place
# (urban, region, cluster travel times, population density). Missing values are kept as NA;
# the gradient-boosting model handles them.
#
# DHS data are licensed to this project: the output stays in data/derived and is not shared.
# Outputs: data/derived/slip_contacts_features.csv (one row per contact, with the calendar month
#          the dose fell due), penta1_by_wealth_north.csv

suppressMessages({library(haven); library(dplyr)})
kr <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA") |>
  mutate(age_m = v008 - b3, region = as.character(as_factor(v024))) |>
  filter(b5 == 1, age_m >= 12, age_m <= 35)

day <- function(d, m, y) {
  ok <- !is.na(d) & d >= 1 & d <= 31 & m >= 1 & m <= 12 & y >= 2015 & y <= 2023
  out <- rep(as.Date(NA), length(d))
  out[ok] <- as.Date(sprintf("%04d-%02d-%02d", as.integer(y[ok]), as.integer(m[ok]), as.integer(d[ok])),
                     optional = TRUE)
  out
}
num <- function(x, bad = numeric(0)) { x <- as.numeric(x); x[x %in% bad] <- NA; x }
kr <- kr |> mutate(
  dob = day(b17, b1, b2),
  p1 = as.Date(ifelse(h3 == 1, day(h3d, h3m, h3y), NA)),
  p2 = as.Date(ifelse(h5 == 1, day(h5d, h5m, h5y), NA)),
  p3 = as.Date(ifelse(h7 == 1, day(h7d, h7m, h7y), NA))) |>
  filter(!is.na(p1), !is.na(dob))

base <- kr |> transmute(
  child_id = paste(caseid, bidx), cluster = v001, region, north = region %in% c("northern", "savannah", "north east", "upper east", "upper west"),
  weight = v005 / 1e6,
  female = as.integer(b4 == 2), birth_order = num(bord), birth_month = as.numeric(b1),
  twin = as.integer(b0 > 0),
  mother_age = num(v012), mother_edu = num(v106, 9), mother_working = as.integer(v714 == 1),
  insured = as.integer(v481 == 1),
  media = rowSums(cbind(num(v157, 9), num(v158, 9), num(v159, 9)) >= 1, na.rm = TRUE),
  problem_permission = as.integer(v467b == 1), problem_money = as.integer(v467c == 1),
  problem_distance = as.integer(v467d == 1), problem_alone = as.integer(v467f == 1),
  anc_visits = num(m14, c(98, 99)), facility_delivery = as.integer(num(m15, 99) >= 20 & num(m15, 99) < 96),
  wealth = num(v190), hh_size = num(v136), under5_in_hh = num(v137),
  head_female = as.integer(v151 == 2), urban = as.integer(v025 == 1),
  g12 = ifelse(!is.na(p2) & p2 > p1, as.numeric(p2 - p1), NA),
  g23 = ifelse(!is.na(p2) & !is.na(p3) & p3 > p2, as.numeric(p3 - p2), NA),
  miss2 = h5 == 0, miss3 = h5 == 0 | (h7 == 0 & !is.na(p2)),
  # Calendar date each contact fell due (schedule 6, 10, 14 weeks: 4 weeks after the previous dose).
  due2 = p1 + 28, due3 = as.Date(ifelse(!is.na(p2), p2 + 28, p1 + 56)))

slip <- function(miss, g) ifelse(miss, 1, ifelse(is.na(g), NA, as.integer(g > 56)))
stype <- function(miss, g) ifelse(miss, "never given", ifelse(is.na(g), NA, ifelse(g > 56, "late", "timely")))
rows <- bind_rows(
  base |> mutate(contact = 2, slip = slip(miss2, g12), slip_type = stype(miss2, g12), due = due2),
  base |> mutate(contact = 3, slip = slip(miss3, g23), slip_type = stype(miss3, g23), due = due3)) |>
  filter(!is.na(slip)) |>
  mutate(due_month = as.integer(format(due, "%m")), due_year = as.integer(format(due, "%Y")),
         smc_season = as.integer(due_month %in% 6:9)) |>
  select(-g12, -g23, -miss2, -miss3, -due2, -due3, -due)

clus <- read.csv("data/derived/dhs_clusters.csv")[, c("cluster", "tt_walk_any", "tt_motor_any", "logpopdens")]
rows <- left_join(rows, clus, by = "cluster")
write.csv(rows, "data/derived/slip_contacts_features.csv", row.names = FALSE)

# Penta1 (any source: card date, card mark or mother's report) by wealth quintile, northern
# children 12-23 months, survey-weighted: the first-contact reach of routine services, for the
# equity analysis (30_equity.py).
p1w <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA") |>
  mutate(age_m = v008 - b3, region = as.character(as_factor(v024)), w = v005 / 1e6) |>
  filter(b5 == 1, age_m >= 12, age_m <= 23,
         region %in% c("northern", "savannah", "north east", "upper east", "upper west")) |>
  group_by(wealth = as.numeric(v190)) |>
  summarise(penta1 = weighted.mean(h3 %in% 1:3, w), n = n(), .groups = "drop")
write.csv(p1w, "data/derived/penta1_by_wealth_north.csv", row.names = FALSE)
print(p1w)
cat("contacts:", nrow(rows), "| children:", nrow(base), "| north contacts:", sum(rows$north),
    "| slip rate all", round(mean(rows$slip), 3), "north", round(mean(rows$slip[rows$north]), 3), "\n")
cat("missing share by feature:\n"); print(round(colMeans(is.na(rows)), 2)[colMeans(is.na(rows)) > 0])
