# WEA Cetacean Baseline Report

R analysis pipeline and report for summarising cetacean baseline information for
the Wind Energy Area (WEA) study region on the Eastern Scotian Shelf.

The project combines vessel sightings from WSDB, DFO aerial survey records, PAM
detections, and spatial boundary layers to produce:

- deduplicated and harmonised cetacean sightings tables
- spatial coverage, confidence, effort, and seasonality summaries
- species and group-level sightings maps
- kernel density estimate (KDE) rasters, contours, and comparison maps
- an HTML report assembled from saved figures and lightweight CSV exports

## Repository Structure

```text
.
├── RUN_PIPELINE.R                 # Master pipeline controller
├── Cetacean_dist_report.Rmd       # RMarkdown report source
├── Cetacean_dist_report.html      # Rendered report
├── scripts/
│   ├── setup/                     # One-time setup and source data cleaning
│   ├── helpers/                   # Shared spatial and grid helpers
│   └── *.R                        # Numbered analysis stages
├── input/
│   ├── raw_data/                  # Raw WSDB, aerial, PAM, and lookup inputs
│   └── processed_data/            # Processed intermediate inputs
├── shapefiles/                    # Spatial boundary and reference layers
└── output/
    ├── data/                      # CSV/RDS summary outputs
    ├── figs/                      # PNG report figures
    ├── shapes/                    # Contour shapefiles
    └── tif/                       # KDE raster outputs
```

## Requirements

This project is written in R and is intended to be run from the project root,
preferably by opening `WEA-cetacean-baseline-report.Rproj` in RStudio.

Core R packages used across the scripts include:

- `tidyverse`, `dplyr`, `tidyr`, `readr`, `stringr`, `lubridate`
- `sf`, `terra`, `spatstat`
- `ggplot2`, `ggspatial`, `patchwork`, `RColorBrewer`, `viridis`
- `here`, `pacman`, `readxl`, `kableExtra`, `knitr`, `rmarkdown`
- `smoothr`, `data.table`, `parallel`, `grafify`, `leaflet`, `scales`

Most scripts use `pacman::p_load()` and will attempt to install missing
packages when run interactively. Spatial packages such as `sf` and `terra` may
also require system geospatial libraries depending on your R installation.

## Key Inputs

Important input paths are configured inside the numbered scripts and spatial
helper files. Common inputs include:

- `input/raw_data/aerial/AllPlatforms_DFO_systematic_surveys_SS_BF_CS_Geoloc-20251205.csv`
- `input/processed_data/aerial/deduplicated_sightings.csv`
- `input/raw_data/wsdb/WSDB_DFO_Sep2025.csv`
- `input/raw_data/PAM/baleen_presence_laura_2025.csv`
- `input/raw_data/speciesCodes.csv`
- `input/raw_data/species_lookup_aerial_to_wsdb.csv`
- spatial layers listed in `scripts/helpers/helper_spatial.R`

To update shapefile locations, edit the `SPATIAL_PATHS` list in
`scripts/helpers/helper_spatial.R`.

## Running the Pipeline

Use `RUN_PIPELINE.R` as the main entry point. The script contains stage toggles
near the top and runs only the stages set to `TRUE`.

From R:

```r
source("RUN_PIPELINE.R")
```

From a terminal:

```sh
Rscript RUN_PIPELINE.R
```

Before running, review the toggles in `RUN_PIPELINE.R`. The controller includes
presets for common workflows such as a full rebuild, a new data delivery, and
map-layout-only updates.

### Pipeline Stages

| Stage | Toggle | Script | Purpose |
| --- | --- | --- | --- |
| Setup | `RUN_SETUP_STUDY_AREA` | `scripts/setup/02_build_study_area.R` | Build the study area polygon from boundary layers |
| Setup | `RUN_CLEAN_AERIAL` | `scripts/setup/01_clean_aerial_sightings.R` | Clean and deduplicate DFO aerial sightings |
| Data prep | `RUN_COMBINE_SIGHTINGS` | `scripts/02_prepare_combined_sightings.R` | Combine WSDB and aerial sightings into one deduplicated dataset |
| Plots | `RUN_PLOT_SIGHTINGS` | `scripts/03_plot_combined_sightings.R` | Create source and species sightings maps |
| Summaries | `RUN_SUMMARIZE_EFFORT` | `scripts/04_summarize_wsdb_effort.R` | Build gridded effort, coverage, and seasonal summaries |
| Summaries | `RUN_SUMMARIZE_CONFIDENCE` | `scripts/05_summarize_wsdb_confidence.R` | Summarise confidence metrics for WEA and study-region cells |
| Summaries | `RUN_SUMMARIZE_SEASONALITY` | `scripts/06_summarize_wsdb_seasonality.R` | Create seasonal cetacean coverage outputs |
| Summaries | `RUN_SUMMARIZE_PAM` | `scripts/07_summarize_pam_detections.R` | Summarise PAM detections by species and season |
| KDE | `RUN_KDE_MODELS` | `scripts/08_run_kde_models.R` | Compute KDE rasters and maps |
| KDE | `RUN_KDE_MAPS_ONLY` | `scripts/09_plot_kde_maps.R` | Recreate KDE map PNGs from existing rasters |
| KDE | `RUN_KDE_CONTOURS` | `scripts/10_create_kde_contours.R` | Generate multi-quantile contour shapefiles from KDE rasters |
| Comparison | `RUN_COMPARE_PAM_SIGHTINGS` | `scripts/11_compare_kde_pam_vs_sightings.R` | Compare PAM and sightings KDE products |

`RUN_KDE_MODELS` can take 20 to 40 minutes. If KDE rasters already exist in
`output/tif/` and only map styling changed, use `RUN_KDE_MAPS_ONLY` instead.
The pipeline stops if `RUN_KDE_MODELS` and `RUN_KDE_MAPS_ONLY` are both enabled.

## Rendering the Report

The report is built from saved PNGs and CSV summaries, so the relevant pipeline
stages should be run before knitting.

```r
rmarkdown::render("Cetacean_dist_report.Rmd")
```

The rendered output is written to `Cetacean_dist_report.html`.

## Outputs

Common generated outputs include:

- `output/data/combined_dedup_1km_day.csv`
- `output/data/deduplication_summary.csv`
- `output/data/sightings_summary.csv`
- `output/data/seasonal_summary_WEA_vs_study.csv`
- `output/data/pam_summary.csv`
- `output/figs/sights/*.png`
- `output/figs/Effort_maps/*.png`
- `output/figs/KDE_maps/*.png`
- `output/figs/compare_sightings_pam/*.png`
- `output/tif/*.tif`
- `output/shapes/multi/*.shp`

## Notes for Future Updates

- Run scripts from the project root so `here::here()` resolves paths correctly.
- Keep `scripts/helpers/helper_spatial.R` as the single source of truth for
  spatial layer paths and CRS settings.
- Run `04_summarize_wsdb_effort.R` before
  `05_summarize_wsdb_confidence.R`; the confidence script depends on objects
  created by the effort script in the same R session.
- Re-run `02_prepare_combined_sightings.R` whenever WSDB, aerial sightings, or
  species lookup inputs change.
- Re-run KDE models only when sightings data, PAM data, bandwidth, resolution,
  or species grouping choices change.
