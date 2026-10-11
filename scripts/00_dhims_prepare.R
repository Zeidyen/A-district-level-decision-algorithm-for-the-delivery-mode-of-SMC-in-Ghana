# DHIMS2 district-month malaria data, all Ghana (Zeid, 9 Oct 2026):
# data/raw/dhims/Malaria_population_district_data.xlsx, sheet "National": 248 districts x 132
# months (2015-2025), same 82 columns as the earlier regional extracts. Checked: the
# Savannah / Northen / North East sheets and the earlier Upper East / Northern / Upper West
# files match the National sheet exactly. The sheet is padded with empty rows (dimension
# to row 1,048,576); those are dropped. Three district-months have zero confirmed cases with
# no tests recorded (early-2015 reporting gaps): Ayawaso Central 2015-11, Yunyoo-Nasuan
# 2015-01, North Gonja 2015-03. They are kept as reported (outside the 2020-24 window used).
# Output: data/raw/dhims/dhims_national_2015_2025.csv
suppressMessages({library(readxl)})
x <- read_excel("data/raw/dhims/Malaria_population_district_data.xlsx", sheet = "National")
x <- x[rowSums(!is.na(x)) > 0, ]
stopifnot(nrow(x) == 248 * 132, length(unique(x$District)) == 248)
x$m <- match(x$month, month.name)
# Key to the geoBoundaries ADM2 names (letters only, lower case), with the spellings that
# differ between DHIMS and geoBoundaries mapped explicitly. Checked for the 55 northern
# districts below.
norm <- function(v) gsub("[^a-z]", "", tolower(v))
alias <- c(gushiegu = "gushegu", tatalesangule = "tatalesanguli", sagnarigu = "sagnerigu",
           tamale = "tamalemetropolitan", yendi = "yendimunicipal", lambussie = "lambussiekarni",
           bolgatangaeast = "bolgaeast", kasenanankana = "kasenanankanaeast",
           westmamprusi = "westmamprusimunicipal", eastgonja = "eastgonjamunicipal")
x$adm2_key <- norm(x$District)
x$adm2_key <- ifelse(x$adm2_key %in% names(alias), alias[x$adm2_key], x$adm2_key)
suppressMessages(library(sf)); sf_use_s2(FALSE)
NORTH <- c("Northern", "Savannah", "North East", "Upper East", "Upper West")
adm1 <- st_read("data/raw/gha_ADM1.geojson", quiet = TRUE)
adm2 <- st_read("data/raw/gha_ADM2.geojson", quiet = TRUE)
adm2$region <- sub(" Region$", "", adm1$shapeName[sapply(st_intersects(st_point_on_surface(adm2), adm1), `[`, 1)])
nk <- norm(adm2$shapeName[adm2$region %in% NORTH])
dk <- unique(x$adm2_key[x$Region %in% NORTH])
stopifnot(length(dk) == 55, setequal(dk, nk))
write.csv(x, "data/raw/dhims/dhims_national_2015_2025.csv", row.names = FALSE)
cat("wrote", nrow(x), "district-months,", length(unique(x$District)), "districts,",
    length(unique(x$Region)), "regions\n")
