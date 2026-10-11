#!/bin/zsh
# WorldPop 100 m covariates for Ghana used by 24_ml_covariates.R (night-time lights, OSM distances,
# slope, ESA CCI land-cover distances). Run from the repository root.
set -e
mkdir -p data/raw/worldpop_cov && cd data/raw/worldpop_cov
B=https://data.worldpop.org/GIS/Covariates/Global_2000_2020/GHA
for p in VIIRS/gha_viirs_100m_2016.tif OSM/DST/gha_osm_dst_road_100m_2016.tif \
         OSM/DST/gha_osm_dst_roadintersec_100m_2016.tif OSM/DST/gha_osm_dst_waterway_100m_2016.tif \
         Slope/gha_srtm_slope_100m.tif ESA_CCI_Annual/2015/gha_esaccilc_dst190_100m_2015.tif \
         ESA_CCI_Annual/2015/gha_esaccilc_dst011_100m_2015.tif ESA_CCI_Annual/2015/gha_esaccilc_dst040_100m_2015.tif; do
  [ -f "$(basename $p)" ] || curl -fsSL -O "$B/$p"
done
