# Download Malaria Atlas Project layers for Ghana into data/raw/map/.
# Friction surfaces (step 1), MAP's own travel time to healthcare (benchmark for
# step 1), and the 2025 Pf prevalence and incidence-rate surfaces (step 3).
suppressMessages({library(malariaAtlas); library(terra); library(sf)})
dir.create("data/raw/map", showWarnings = FALSE, recursive = TRUE)
gha <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
ext <- st_as_sf(st_as_sfc(st_bbox(gha)))

get <- function(id, out, year = NULL) {
  f <- file.path("data/raw/map", out)
  if (file.exists(f)) return(invisible(cat("have", out, "\n")))
  r <- getRaster(dataset_id = id, shp = ext, year = year)
  writeRaster(r, f, overwrite = TRUE)
  cat("got", out, "\n")
}

get("Accessibility__202001_Global_Walking_Only_Friction_Surface", "friction_walking_2019.tif")
get("Accessibility__202001_Global_Motorized_Friction_Surface", "friction_motorised_2019.tif")
get("Accessibility__202001_Global_Walking_Only_Travel_Time_To_Healthcare", "map_tt_walking_healthcare.tif")
get("Accessibility__202001_Global_Motorized_Travel_Time_to_Healthcare", "map_tt_motorised_healthcare.tif")
get("Malaria__202608_Global_Pf_Parasite_Rate", "pfpr_2025.tif", year = 2025)
get("Malaria__202608_Global_Pf_Incidence_Rate", "pf_incidence_rate_2025.tif", year = 2025)
