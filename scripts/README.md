```markdown
# Analysis scripts

This directory contains the R scripts used for data processing, statistical analysis, extrapolation, sensitivity analyses, and generation of the main and supplementary figures and tables.

The scripts are numbered according to the approximate order of the analysis workflow.

## Workflow overview

### 1. Final_00–Final_09: Upstream data processing and model preparation

These scripts document the upstream processing used to generate the analysis-ready datasets and checkpoint files.

The upstream workflow uses data from multiple sources, including:

- Demographic and Health Surveys (DHS)
- Institute for Health Metrics and Evaluation (IHME) Global Burden of Disease (GBD) data
- Antibiotic consumption data
- World Bank and other publicly available country-level covariates

Some of the source data required for these scripts cannot be redistributed through this repository. In particular, individual-level DHS microdata and specific IHME/GBD source-data extracts are not included.

Therefore, scripts `Final_00`–`Final_09` are provided primarily to document the upstream analytical workflow. Running the complete upstream workflow requires users to independently obtain the required source data and reconstruct the necessary source-data extracts.

See `data/README.md` and `data/datasets/README.md` for additional information.

### 2. Final_10–Final_19: Publicly reproducible downstream main analyses

The publicly supported reproduction workflow for the main analysis begins with `Final_10`.

The required analysis checkpoint files are archived separately on Zenodo. After downloading these files and placing them in the repository's `outputs/` directory, the downstream main analyses can be run beginning with:

`Final_10_Extrapolation_and_absolute_number_revised.R`

The subsequent scripts generate the downstream estimates, tables, and figures used in the manuscript and supplementary materials.

See `outputs/README.md` for the list of required checkpoint files and instructions.

### 3. Final_20: Full sensitivity-analysis workflow

`Final_20_Sensitivity_analysis.R` documents the complete sensitivity-analysis workflow, including the upstream processing used to generate sensitivity-specific analysis inputs and checkpoints.

The sensitivity analysis begins from individual-level DHS-derived data and also uses specific IHME/GBD source-data extracts. Because these upstream inputs cannot be redistributed through this repository, `Final_20` is provided primarily for transparency and documentation of the complete sensitivity-analysis workflow.

Running the entire script requires users to independently obtain and reconstruct the required upstream source-data inputs.

### 4. Final_21: Publicly reproducible sensitivity analysis

The publicly supported sensitivity-analysis workflow begins with `Final_21`.

Sensitivity-specific checkpoint files are archived on Zenodo together with the main-analysis checkpoint files. After downloading the checkpoint files and placing them in the appropriate locations under `outputs/`, the downstream sensitivity analyses can be run beginning with:

`Final_21_Sensitivity_analysis_from_checkpoint.R`

This script reproduces the downstream sensitivity-analysis estimates, tables, and figures without requiring access to the individual-level DHS microdata or original IHME/GBD source-data extracts.

See `outputs/README.md` for the required sensitivity-analysis checkpoint files and instructions.

## Running the scripts

The repository root should be used as the working directory so that relative paths such as `data/`, `outputs/`, and `figures/` resolve correctly.

### Main analysis

For the public reproduction of the main analysis:

1. Download the Zenodo checkpoint files.
2. Place them in the locations specified in `outputs/README.md`.
3. Run the scripts beginning with `Final_10_Extrapolation_and_absolute_number_revised.R`.
4. Continue through `Final_19` as needed to reproduce the downstream analyses and figures.

### Sensitivity analysis

For the public reproduction of the sensitivity analysis:

1. Download the sensitivity-analysis checkpoint files from the same Zenodo archive.
2. Place them in the locations specified in `outputs/README.md`.
3. Run `Final_21_Sensitivity_analysis_from_checkpoint.R`.

`Final_20_Sensitivity_analysis.R` should only be run when the required upstream DHS and IHME/GBD source-data inputs are available.

Scripts may create additional files and subdirectories within `outputs/`.

## Software

The analyses were conducted in R. Required R packages are loaded or installed within the individual analysis scripts.

Additional software and package-version information will be provided in the main repository README.
```