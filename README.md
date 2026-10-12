# A district-level decision algorithm for the delivery mode of SMC in Ghana

Code for the analysis *"A district-level decision algorithm for the delivery mode of seasonal malaria chemoprevention: development and application in northern Ghana"* (Anafo AY et al., manuscript in preparation).

The analysis asks where routine delivery of seasonal malaria chemoprevention (SMC), through health facilities and community-based health planning and services (CHPS) compounds, could protect as many children as the door-to-door campaign in the five SMC regions of Ghana (55 districts). It combines access, a Bayesian spatial model of late or missed routine contacts, observed malaria burden, and coverage under each delivery mode into a six-step district decision algorithm:

1. **Campaign check** — is survey-measured campaign coverage below the programme target?
2. **Age reach** — does the routine model reach children aged 2–4 years?
3. **Reliability** — break-even k (how much more often a routine SMC visit could be missed than a vaccination visit) against the evidence for k, with its uncertainty.
4. **Equity safeguard** — would the poorest wealth quintile bear significantly more than its share of the protection lost?
5. **Pilot priority** — ranking of pilot districts.
6. **Measure and decide** — what a pilot should measure, and the rule for switching or reverting.

**No data are included in this repository.** See [Data](#data) for how to obtain each input.

## Running the pipeline

```bash
Rscript install_packages.R           # R packages (INLA from its own repository)
pip install -r requirements.txt      # Python packages
zsh scripts/download_worldpop_covariates.sh
./run_all.sh                         # whole pipeline; logs in out/logs/
```

Run everything from the repository root. Inputs go in `data/raw/` (layout in [`data/raw/README.md`](data/raw/README.md)); intermediate results are written to `data/derived/` and figures to `out/`. Posterior sampling is seeded and single-threaded, so results are reproducible.

## Scripts

| Script | Purpose |
| --- | --- |
| `00_dhims_prepare.R` | DHIMS2 district-month malaria data: clean and match districts to boundaries |
| `01_download_map.R` | Download Malaria Atlas Project friction surfaces |
| `01b_facilities.R` | Clean the facility list (public and faith-based facilities; CHPS flag) |
| `02_access.R`, `05_access_ghana.R` | Least-cost travel time to the nearest facility and CHPS |
| `03_dhs_clusters.R`, `04_timeliness.R` | DHS 2022 cluster covariates (population-weighted buffers) and per-contact slip outcomes |
| `06_inla_timeliness.R`, `lib_inla.R` | INLA-SPDE model of late or missed routine contacts (step 2) |
| `07a_chirps.R`, `07_rain_timeliness.R` | CHIRPS rainfall and its association with slips |
| `08_checks.R` | Model decomposition, spatial block cross-validation |
| `09_burden.R` | Under-5 malaria burden in the months SMC protects (step 3) |
| `10_step4_impact.R` | Coverage and cases averted under each delivery mode; model variants (step 4) |
| `11_summary_map.R`, `12_k_chart.R`, `14_sensitivity.R` | Robustness across model variants |
| `15_wet_access.R`, `16_env_plots.R`, `17_terrain.R` | Rainy-season flooding, environmental and terrain checks |
| `18_age_reach.R` | Age-reach factor from vitamin A receipt by age |
| `19_dhims_qa.R`, `20_validate_fingerprint.R` | DHIMS2 quality checks and before/after SMC check |
| `21_breakeven_gsa.R` | Break-even k, break-even campaign coverage, global sensitivity analysis |
| `22_district_cv.R` | Leave-one-district-out cross-validation of the spatial model |
| `24_ml_covariates.R`, `25_ml_stack.py`, `26_ml_stack_cv.R` | Stacked machine-learning covariate model and its cross-validation |
| `27a_slip_features.R`, `27_slip_drivers.py` | Predictors of slipped contacts (gradient boosting, SHAP) |
| `28_season_slips.R` | Do contacts slip more in the SMC season? |
| `29_hybrid.R` | Hybrid delivery: children moved against cases averted lost |
| `30_equity.py` | Who would lose protection (wealth quintiles, regional disproportionality ratios) |
| `31_delivery_algorithm.R` | The SMC delivery-mode algorithm, its map and flowchart |
| `32_k_dropout_bias.R` | Bias from calibrating k on dropout rather than lateness |
| `33_sobol_total.R` | First-order and total-effect Sobol indices |
| `34_study_site_map.R` | Study site map: the five SMC regions and their 55 districts, with an Africa inset |
| `13_summary_docx.py`, `23_manuscript_docx.py` | Report generators (results summary; Word version of the manuscript) |
| `compare_runs.py` | Compare outputs between two runs |

The algorithm's thresholds (evidence range for k, confidence level, campaign target, equity rule) are settings at the top of `31_delivery_algorithm.R` and can be changed by a programme.

## Data

| Input | Source | Access |
| --- | --- | --- |
| Ghana DHS 2022 (child recode KR, GPS GE) | [DHS Program](https://dhsprogram.com/) | Free, on registration and approval of a research request; may not be redistributed |
| DHIMS2 district-month malaria data, 2015–2025 | Ghana Health Service | On request to the Ghana Health Service; not public |
| Public health facility list with coordinates | [Ghana Health Service facility finder](https://ghs.gov.gh/facility-finder) | Public website; no stated licence, so not redistributed here |
| Maina et al. public facility list (sensitivity variant) | [HDX](https://data.humdata.org/dataset/health-facilities-in-sub-saharan-africa) | Open |
| Friction surfaces | [Malaria Atlas Project](https://malariaatlas.org/) | Open; downloaded by `01_download_map.R` |
| Population by age, R2025A, 1 km | [WorldPop](https://hub.worldpop.org/) | Open (CC BY 4.0) |
| Satellite and access covariates | WorldPop covariates | Open; `scripts/download_worldpop_covariates.sh` |
| Monthly rainfall | [CHIRPS v2.0](https://data.chc.ucsb.edu/products/CHIRPS-2.0/) | Open; `07a_chirps.R` |
| Surface water seasonality | [JRC Global Surface Water](https://global-surface-water.appspot.com/) | Open |
| Elevation | SRTM via the `geodata` R package | Open |
| Administrative boundaries | [geoBoundaries](https://www.geoboundaries.org/) gbOpen GHA ADM1 and ADM2 | Open |
| Country outlines (map inset) | GADM, via the `geodata` R package | Free for academic use; downloaded by `34_study_site_map.R` |
| SMC administrative coverage (optional) | National Malaria Elimination Programme | On request; used only if supplied |

## Citation

Please cite the paper once published. Citation details will be added here.

## Licence

The code is released under the [MIT License](LICENSE). The licence covers the code only; the input data keep their own terms (see [Data](#data)).

## Contact

Abdulzeid Yen Anafo — University of Mines and Technology (UMaT), Tarkwa, Ghana.
