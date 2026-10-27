```markdown
# Data

This directory contains the publicly shareable input data used in the analysis.

Some source datasets used in the upstream analysis cannot be redistributed through this repository. In particular, individual-level Demographic and Health Surveys (DHS) microdata and selected Institute for Health Metrics and Evaluation (IHME) Global Burden of Disease (GBD) inputs must be obtained directly from the original data providers.

## Included data

The following publicly shareable inputs are included:

- `DDD_country_year_class.csv` — country-year antibiotic consumption data used in the antibiotic consumption analyses.
- `GDP_PPP.csv` — gross domestic product per capita (purchasing power parity) data.
- `Income_level.csv` — country income-group classifications.
- `gbd_region.csv` — country-to-GBD-region classification.
- `population_under_5.csv` — population estimates for children younger than 5 years.
- `Precipitation.csv` — annual precipitation data.
- `Temperature_2024.csv` — temperature data.
- `WHO2024_malarial_endemic.csv` — malaria-endemicity classification.
- `ver.2_res.fev.prop.csv` — literature-review data used to estimate the proportion of respiratory infection among children with fever.
- `UNIndex/` — World Development Indicators used to construct the principal-component socioeconomic index.

Detailed source and citation information for these datasets will be provided in the main repository documentation.

## Data not redistributed

### 1. DHS microdata

Individual-level DHS microdata are not included in this repository.

Users wishing to reproduce the DHS data-processing steps should obtain authorized access directly from The DHS Program and place the prepared survey files in:

- `data/datasets/`

See `data/datasets/README.md` for instructions.

### 2. IHME/GBD inputs

The following upstream IHME/GBD input files used in the analysis are not redistributed here:

- `All_inf2.csv`
- `malaria_proportion.csv`

These files represent specific data extracts used in the original analysis. Reproducing the upstream IHME/GBD processing therefore requires access to the corresponding source data and reconstruction of matching extracts with the required variables and structure.

The scripts documenting these upstream processing steps are provided in the repository, but the upstream analysis is not intended to be fully self-contained without the original source-data extracts.

## Reproduction paths

Two analysis routes are documented.

### 1. Upstream processing from source data

Scripts beginning with:

- `scripts/Final_00_DHS_data.R`

document the processing used to generate the analysis-ready datasets and model inputs from the original DHS, IHME/GBD, and other source data.

Running these upstream scripts requires users to independently obtain the required source datasets and, where applicable, reproduce the source-data extracts used in the original analysis.

### 2. Public reproduction from analysis checkpoints

The publicly supported reproduction workflow begins from analysis checkpoint files archived in the associated Zenodo record.

This workflow does not require redistribution of individual-level DHS microdata or the original IHME/GBD source extracts.

After downloading the checkpoint files into the repository's `outputs/` directory, run the analysis beginning with:

- `scripts/Final_10_Extrapolation_and_absolute_number_revised.R`

This reproduces the downstream analyses, tables, and figures from the publicly archived analysis checkpoints.

See `outputs/README.md` for the required checkpoint files and instructions.
```