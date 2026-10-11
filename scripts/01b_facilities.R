# Facilities for the access step: one cleaned file, whichever source is used.
#
# Source (env var FACILITY_SOURCE):
#   ghs   (default from 10 Oct 2026) Ghana Health Service public facility finder,
#         https://ghs.gov.gh/facility-finder, pulled 2026-10-10 from its open JSON endpoint
#         (data/raw/facilities/ghs_api/; see README there). 10,339 records, 6,974 CHPS
#         nationally, 1,675 in the five northern regions. No licence stated: ask GHS before
#         publishing.
#   maina Maina et al. 2019 public-sector list (HDX), the earlier source: ~632 CHPS nationally.
#
# Kept: government, faith-based (CHAG and other), quasi-government and public facilities;
# routine SMC would run through the public system, so private facilities (and mines) are
# dropped. Records with no coordinates (0,0), outside Ghana, or exact duplicates (same name
# and coordinates) are dropped. Region is assigned from the coordinates (the GHS file has no
# district field; a few points lie just across the region named in the record).
#
# Output: data/derived/facilities.csv (source, name, type, is_chps, ownership, region_stated,
#         region, lon, lat)

suppressMessages({library(sf); library(readxl)})
sf_use_s2(FALSE)
SRC <- Sys.getenv("FACILITY_SOURCE", "ghs")
KEEP_OWN <- c("GOVERNMENT", "CHAG", "OTHER FAITH-BASED", "QUASI-GOVERNMENT", "PUBLIC")

if (SRC == "ghs") {
  x <- read.csv("data/raw/facilities/ghs_api/ghs_public_facilities_20261010.csv")
  f <- data.frame(source = "GHS facility finder 2026-10-10", name = x$name, type = x$facilityType,
                  is_chps = x$facilityType == "CHPS", ownership = x$ownership,
                  region_stated = x$region, lon = x$longitude, lat = x$latitude)
  n0 <- nrow(f)
  f <- f[f$ownership %in% KEEP_OWN, ]; n_own <- nrow(f)
} else if (SRC == "maina") {
  x <- read_excel("data/raw/facilities/ssa_public_facilities.xlsx")
  x <- x[x$Country == "Ghana", ]
  f <- data.frame(source = "Maina et al. 2019 (HDX)", name = x$Facility_n, type = x$Facility_t,
                  is_chps = grepl("Community-based", x$Facility_t), ownership = x$Ownership,
                  region_stated = x$Admin1, lon = x$Long, lat = x$Lat)
  n0 <- nrow(f); n_own <- nrow(f)                        # already public sector (MoH + FBO)
} else stop("FACILITY_SOURCE must be ghs or maina")

f <- f[!(f$lat == 0 | f$lon == 0) & !is.na(f$lat) & !is.na(f$lon), ]; n_xy <- nrow(f)
f <- f[!duplicated(f[, c("name", "lat", "lon")]), ]; n_dup <- nrow(f)
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
p <- st_as_sf(f, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
f$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(p, adm1), function(i) if (length(i)) i[1] else NA)])
f <- f[!is.na(f$region), ]
write.csv(f, "data/derived/facilities.csv", row.names = FALSE)

NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
cat(sprintf("source %s: %d records -> %d public/faith -> %d with coordinates -> %d after de-duplication -> %d inside Ghana\n",
            SRC, n0, n_own, n_xy, n_dup, nrow(f)))
cat("Five northern regions:", sum(f$region %in% NORTH), "facilities,", sum(f$is_chps & f$region %in% NORTH), "CHPS\n")
print(table(f$region[f$region %in% NORTH], f$is_chps[f$region %in% NORTH], dnn = c("region", "CHPS")))
