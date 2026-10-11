# Does routine-contact TIMELINESS track travel time, where coverage did not?
#
# Ghana's schedule gives Penta at 6, 10 and 14 weeks. SMC needs a contact every
# month, so the measure closest to it is whether a scheduled monthly contact
# slipped by four weeks or more. Using card-dated doses, children 12-35 months:
#
#   late_p3   Penta3 given after 18 weeks of age (due at 14, plus 4 weeks' grace)
#   gap_slip  either gap (dose 1->2 or 2->3) longer than 8 weeks, i.e. a monthly
#             contact missed by at least a month
#   gap_days  total days from Penta1 to Penta3 (scheduled: 56)
#
# Children without card dates are excluded. Card possession may itself depend
# on distance, so that is reported as a check.

suppressMessages({library(haven); library(dplyr); library(glmmTMB)})

north <- c("northern", "savannah", "north east", "upper east", "upper west")
kr <- read_dta("data/raw/dhs/GHKR8CDT/GHKR8CFL.DTA") |>
  mutate(age_m = v008 - b3, region = as.character(as_factor(v024))) |>
  filter(b5 == 1, age_m >= 12, age_m <= 35)    # all Ghana; north selected below

day <- function(d, m, y) {
  ok <- !is.na(d) & d >= 1 & d <= 31 & m >= 1 & m <= 12 & y >= 2015 & y <= 2023
  out <- rep(as.Date(NA), length(d))
  out[ok] <- as.Date(sprintf("%04d-%02d-%02d", as.integer(y[ok]), as.integer(m[ok]), as.integer(d[ok])),
                     optional = TRUE)
  out
}
kr <- kr |> mutate(
  dob = day(b17, b1, b2),
  p1 = ifelse(h3 == 1, day(h3d, h3m, h3y), NA) |> as.Date(),
  p2 = ifelse(h5 == 1, day(h5d, h5m, h5y), NA) |> as.Date(),
  p3 = ifelse(h7 == 1, day(h7d, h7m, h7y), NA) |> as.Date(),
  has_card_dates = !is.na(p1),
  age_p3 = as.numeric(p3 - dob),
  g12 = as.numeric(p2 - p1), g23 = as.numeric(p3 - p2),
  ok = !is.na(dob) & !is.na(p1) & !is.na(p2) & !is.na(p3) &
       age_p3 >= 70 & g12 > 0 & g23 > 0,                 # drop impossible dates
  late_p3  = as.integer(age_p3 > 126),
  gap_slip = as.integer(g12 > 56 | g23 > 56),
  gap_days = g12 + g23)

clus <- read.csv("data/derived/dhs_clusters.csv")[, c("cluster", "urban", "tt_walk_any", "tt_motor_any")]
kr <- left_join(kr, clus, by = c("v001" = "cluster")) |>
  mutate(log2_tt = log2(pmax(tt_walk_any, 5)), wealth = as.integer(v190))

cat("Northern children 12-35 m:", sum(kr$region %in% north),
    "| with card date for Penta1:", sum(kr$has_card_dates[kr$region %in% north]),
    "| all three Penta dates usable:", sum(kr$ok[kr$region %in% north]), "\n\n")

d_all <- filter(kr, ok)                         # all Ghana, for the step 2 model
kr_all <- kr
kr <- filter(kr, region %in% north)
d <- filter(d_all, region %in% north)
cat("(north only from here; Ghana-wide usable children:", nrow(d_all), ")\n")
cat(sprintf("Penta3 after 18 weeks: %.0f%%   any gap > 8 weeks: %.0f%%   median Penta1->3: %.0f days (scheduled 56)\n\n",
            100 * mean(d$late_p3), 100 * mean(d$gap_slip), median(d$gap_days)))

# By walking-time band (child level, unweighted).
d$band <- cut(d$tt_walk_any, c(0, 30, 60, 120, Inf), c("<30", "30-60", "60-120", ">120"), right = FALSE)
print(d |> group_by(band) |> summarise(n = n(),
  late_p3 = round(100 * mean(late_p3)), gap_slip = round(100 * mean(gap_slip)),
  median_gap = median(gap_days)) |> as.data.frame(), row.names = FALSE)

# Mixed models: cluster random intercept, region fixed effects.
# Odds ratio per doubling of walking time, crude and adjusted for wealth + urban.
or <- function(f) {
  # glmmTMB rather than lme4: glmer failed to converge on the crude late_p3 model.
  m <- glmmTMB(f, data = d, family = binomial)
  b <- fixef(m)$cond["log2_tt"]; se <- sqrt(vcov(m)$cond["log2_tt", "log2_tt"])
  sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
}
cat("\nOdds ratio per doubling of walking time [crude | adj. wealth + urban]:\n")
for (y in c("late_p3", "gap_slip")) {
  cat(sprintf("  %-9s %s | %s\n", y,
      or(as.formula(paste(y, "~ log2_tt + region + (1 | v001)"))),
      or(as.formula(paste(y, "~ log2_tt + wealth + urban + region + (1 | v001)")))))
}

# Selection check: are card dates rarer far from facilities?
cat("\nShare of children with Penta1 card date, by band:\n")
kr$band <- cut(kr$tt_walk_any, c(0, 30, 60, 120, Inf), c("<30", "30-60", "60-120", ">120"), right = FALSE)
print(kr |> group_by(band) |> summarise(n = n(), card = round(100 * mean(has_card_dates))) |>
        as.data.frame(), row.names = FALSE)

# ---- contact-level outcome for the step 2 model ----------------------------------
# Each scheduled monthly contact after Penta1 (dose 2, dose 3) is one trial, so the
# model estimates the per-contact slip probability q directly. (The child-level
# "either gap" outcome cannot be converted to q by 1 - sqrt(1 - s): the two gaps
# are correlated within a child, within-child OR 3.2, so that formula underestimated
# q, 0.131 vs 0.144 per gap; found by the 1 Oct audit.)
# Children with a card-dated Penta1 are included even if they never returned: a
# dose never given (h5 / h7 == 0) is a missed contact, i.e. slipped. A dose reported
# only by recall (no card date) leaves that contact unknown (NA, not counted).
# Thresholds: slip = gap > 56 days (main: a monthly contact missed by >= 4 weeks);
# > 49 and > 42 days as sensitivity (partial loss of a 28-day protection window).
contacts <- kr_all |>
  filter(has_card_dates, !is.na(dob)) |>
  transmute(v001,
            g12 = ifelse(!is.na(p2) & p2 > p1, as.numeric(p2 - p1), NA),
            g23 = ifelse(!is.na(p2) & !is.na(p3) & p3 > p2, as.numeric(p3 - p2), NA),
            miss2 = h5 == 0, miss3 = h5 == 0 | (h7 == 0 & !is.na(p2)))
slip_at <- function(thr) {
  c2 <- ifelse(contacts$miss2, 1, ifelse(is.na(contacts$g12), NA, as.integer(contacts$g12 > thr)))
  c3 <- ifelse(contacts$miss3, 1, ifelse(is.na(contacts$g23), NA, as.integer(contacts$g23 > thr)))
  cbind(c2, c3)
}
cc <- data.frame(cluster = contacts$v001)
m56 <- slip_at(56); m49 <- slip_at(49); m42 <- slip_at(42)
cc$contacts <- rowSums(!is.na(m56))
cc$slip56 <- rowSums(m56, na.rm = TRUE); cc$slip49 <- rowSums(m49, na.rm = TRUE)
cc$slip42 <- rowSums(m42, na.rm = TRUE)

# Cluster-level counts for the step 2 model (child-level columns kept for reference).
cl <- d_all |> group_by(cluster = v001) |>
  summarise(n_timely = n(), late_p3 = sum(late_p3), gap_slip = sum(gap_slip), .groups = "drop") |>
  full_join(cc |> group_by(cluster) |>
              summarise(contacts = sum(contacts), slip56 = sum(slip56), slip49 = sum(slip49),
                        slip42 = sum(slip42), .groups = "drop"), by = "cluster")
cl[is.na(cl)] <- 0
write.csv(cl, "data/derived/dhs_clusters_timeliness.csv", row.names = FALSE)

nk <- cl$cluster %in% unique(kr$v001)
cat(sprintf("\nContact-level (all Ghana): %d contacts in %d clusters; per-contact slip >56 d %.3f\n",
            sum(cl$contacts), sum(cl$contacts > 0), sum(cl$slip56) / sum(cl$contacts)))
cat(sprintf("North: %d contacts; slip >56 d %.3f | >49 d %.3f | >42 d %.3f\n", sum(cl$contacts[nk]),
            sum(cl$slip56[nk]) / sum(cl$contacts[nk]), sum(cl$slip49[nk]) / sum(cl$contacts[nk]),
            sum(cl$slip42[nk]) / sum(cl$contacts[nk])))
