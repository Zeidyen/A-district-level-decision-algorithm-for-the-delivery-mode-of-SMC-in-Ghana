# Input data layout

No data are included in this repository. Place the inputs here with these names before running `run_all.sh`:

```
data/raw/
  gha_ADM1.geojson, gha_ADM2.geojson                     geoBoundaries gbOpen, Ghana regions and districts
  gha_t_00_2025_CN_1km_R2025A_UA_v1.tif                 WorldPop R2025A, age <1, 1 km
  gha_t_01_2025_CN_1km_R2025A_UA_v1.tif                 WorldPop R2025A, ages 1-4, 1 km
  dhs/GHKR8CDT/GHKR8CFL.DTA                             Ghana DHS 2022 children's recode (Stata)
  dhs/GHGE8AFL/GHGE8AFL.shp                             Ghana DHS 2022 cluster GPS
  dhims/Malaria_population_district_data.xlsx           DHIMS2 district-month malaria data (sheet "National")
  facilities/ghs_api/ghs_public_facilities_<date>.csv   GHS facility finder export (see 01b_facilities.R for columns)
  facilities/ssa_public_facilities.xlsx                 Maina et al. list (sensitivity variant only)
  map/                                                  filled by 01_download_map.R
  chirps/                                               filled by 07a_chirps.R
  jrc_gsw/                                              JRC Global Surface Water seasonality tiles covering Ghana
  elevation/                                            filled by the geodata package
  worldpop_cov/                                         filled by scripts/download_worldpop_covariates.sh
  nmep/smc_admin_coverage.csv                           optional: region, district, coverage
```

DHS data must be requested from the DHS Program and may not be redistributed. DHIMS2 data are held by the Ghana Health Service.
