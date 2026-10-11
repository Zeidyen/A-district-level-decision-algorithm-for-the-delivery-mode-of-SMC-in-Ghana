# Collate step 4 sensitivity scenarios (run 10_step4_impact.R <scenario> first).
# Writes data/derived/step4_sensitivity.csv: district class counts, confident
# classes, children in red districts and routine shortfall per scenario.
SC <- c(main = "", thr49 = "_thr49", thr42 = "_thr42", walking = "_walking", region = "_region", reach90 = "_reach90",
        age = "_age", adh = "_adh", adh0 = "_adh0", age_adh = "_age_adh",
        nospatial = "_nospatial", regiononly = "_regiononly", stacked = "_stacked")
# Facility-list sensitivity: the main analysis run on the earlier Maina et al. list, kept from
# the last run before the switch to the GHS facility finder (10 Oct 2026), if present.
if (file.exists("data/derived/district_step4_facmaina.csv")) SC <- c(SC, facmaina = "_facmaina")
LAB <- c(main = "Main (gap > 8 wk, motorised, no region effects)", thr49 = "Slip = gap > 7 weeks",
         thr42 = "Slip = gap > 6 weeks", walking = "Walking instead of motorised time",
         region = "Region effects in step 2", reach90 = "Routine reach 90% of Penta1",
         age = "Routine contact falls at 24-59 m (vitamin A)", adh = "Lower routine adherence, partial course 50%",
         adh0 = "Lower routine adherence, partial course 0%", age_adh = "Age + adherence (50%)",
         nospatial = "Step 2 without spatial field", regiononly = "Step 2 = regional means only",
         stacked = "Step 2 = stacked machine-learning covariate + field",
         facmaina = "Earlier facility list (Maina et al., ~212 northern CHPS)")
out <- do.call(rbind, lapply(names(SC), function(s) {
  d <- read.csv(sprintf("data/derived/district_step4%s.csv", SC[s]))
  data.frame(scenario = s, label = LAB[s],
             green = sum(d$class == "green"), amber = sum(d$class == "amber"), red = sum(d$class == "red"),
             green_p80 = sum(d$p_green >= 0.8), red_p80 = sum(d$p_red >= 0.8),
             pct_children_red = round(100 * sum(d$children_3_59m[d$class == "red"]) / sum(d$children_3_59m)),
             shortfall_pct = round(100 * sum(d$shortfall) / sum(d$averted_campaign)))
}))
write.csv(out, "data/derived/step4_sensitivity.csv", row.names = FALSE)
print(out[, -2], row.names = FALSE)

# Robustness of each district's class across all variants -> district_robustness.csv
suppressMessages(library(dplyr))
d <- bind_rows(lapply(names(SC), function(s) read.csv(sprintf("data/derived/district_step4%s.csv", SC[s])) |>
                        mutate(scen = s)))
r <- d |> group_by(region, district) |>
  summarise(n = n(), red = sum(class == "red"), green = sum(class == "green"), amber = sum(class == "amber"),
            ratio_main = ratio[scen == "main"], ratio_min = min(ratio), ratio_max = max(ratio),
            kids = children_3_59m[scen == "main"], .groups = "drop") |>
  mutate(robust = case_when(red == n ~ "red in all", red >= 9 ~ "red in >=9", green == 0 ~ "never green", TRUE ~ "mixed")) |>
  arrange(ratio_main)
write.csv(r, "data/derived/district_robustness.csv", row.names = FALSE)
cat("\nRobustness across", r$n[1], "variants:", paste(names(table(r$robust)), table(r$robust), collapse = " | "),
    "| never green:", sum(r$green == 0), "districts,", round(100 * sum(r$kids[r$green == 0]) / sum(r$kids)), "% of children\n")
