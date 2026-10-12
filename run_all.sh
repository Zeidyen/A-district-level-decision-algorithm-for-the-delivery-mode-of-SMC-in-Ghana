#!/bin/zsh
# Rerun the whole geospatial pipeline in order. Downloads (01, CHIRPS raw) are assumed present.
# Usage: ./run_all.sh  (from geo/). Each step's output goes to out/logs/<step>.log
set -e
mkdir -p out/logs
step() { echo "== $1"; Rscript "scripts/$1" > "out/logs/${1%.*}.log" 2>&1 || { echo "FAILED: $1 (see out/logs/${1%.*}.log)"; exit 1; } }
step 00_dhims_prepare.R    # DHIMS2 national district-month data -> CSV
step 01b_facilities.R      # facilities: GHS facility finder (FACILITY_SOURCE=maina for the old list)
step 02_access.R           # step 1, northern travel time + district tables
step 05_access_ghana.R     # Ghana-wide travel time + density (needed by 03)
step 03_dhs_clusters.R     # DHS cluster covariates (population-weighted buffers)
step 04_timeliness.R       # card-date timeliness, cluster counts
step 06_inla_timeliness.R  # step 2 INLA-SPDE map
step 07a_chirps.R          # CHIRPS stack
step 07_rain_timeliness.R  # rainfall test
step 08_checks.R           # decomposition, spatial CV, sensitivity
step 09_burden.R           # step 3 burden (observed under-5 cases, DHIMS2)
step 18_age_reach.R        # vitamin A by age (input to the "age" scenario)
step 24_ml_covariates.R    # satellite/access covariates for the ML stack (WorldPop 2015-16, CHIRPS, facilities)
echo "== 25_ml_stack.py"; python3 -u scripts/25_ml_stack.py > out/logs/25_ml_stack.log 2>&1 || { echo "FAILED: 25"; exit 1; }
step 10_step4_impact.R     # step 4/5 classification + k stress test (main)
for s in thr49 thr42 walking region reach90 age adh adh0 age_adh nospatial regiononly stacked; do
  echo "== 10_step4_impact.R $s"
  Rscript scripts/10_step4_impact.R $s > out/logs/10_step4_$s.log 2>&1 || { echo "FAILED: scenario $s"; exit 1; }
done
step 14_sensitivity.R      # scenario table + district robustness
step 21_breakeven_gsa.R    # district break-even k + global sensitivity
step 22_district_cv.R      # leave-one-district-out CV: spatial model vs regional means
step 26_ml_stack_cv.R      # stacked ML covariate in INLA; leave-one-district-out vs spatial models
step 27a_slip_features.R   # contact-level DHS features
echo "== 27_slip_drivers.py"; python3 -u scripts/27_slip_drivers.py > out/logs/27_slip_drivers.log 2>&1 || { echo "FAILED: 27"; exit 1; }
step 28_season_slips.R     # do routine contacts slip more in the SMC season (Jun-Sep)?
step 29_hybrid.R           # hybrid delivery: children moved vs cases averted lost
echo "== 30_equity.py"; python3 -u scripts/30_equity.py > out/logs/30_equity.log 2>&1 || { echo "FAILED: 30"; exit 1; }
step 31_delivery_algorithm.R  # SMC delivery-mode algorithm (district decision rule)
step 32_k_dropout_bias.R   # bias from calibrating k on dropout (Supplementary Text S1)
step 33_sobol_total.R      # first-order and total-effect Sobol indices (Supplementary Text S4)
step 34_study_site_map.R   # study site map (greyscale, Africa inset)
echo "== 10_step4_impact.R prior15"; Rscript scripts/10_step4_impact.R prior15 > out/logs/10_step4_prior15.log 2>&1 || { echo "FAILED: prior15"; exit 1; }
step 15_wet_access.R       # rainy-season access (JRC surface water)
step 16_env_plots.R        # environmental layers (rainfall, surface water, elevation, MAP)
step 17_terrain.R          # terrain maps + terrain vs timeliness check
step 19_dhims_qa.R         # DHIMS2 QA + SMC fingerprint
step 20_validate_fingerprint.R
step 11_summary_map.R
step 12_k_chart.R
echo "== 13_summary_docx.py"; python3 scripts/13_summary_docx.py > out/logs/13_summary_docx.log 2>&1
echo "all steps done"
