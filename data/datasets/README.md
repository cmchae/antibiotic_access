```markdown
# DHS microdata

Individual-level Demographic and Health Surveys (DHS) microdata are not redistributed in this repository.

To reproduce the DHS data-processing steps:

1. Register with The DHS Program.
2. Request access to the DHS surveys required for this study.
3. Download the corresponding child recode (KR) datasets.
4. Prepare the survey files in RDS format using the original DHS file names.
5. Place the resulting `.rds` files in this directory (`data/datasets/`).
6. Run `scripts/Final_00_DHS_data.R` to begin the upstream DHS data-processing workflow.

The analysis script automatically reads the `.rds` files placed in this directory.

Subsequent upstream analyses also require additional source data, including specific IHME/GBD data extracts that are not redistributed in this repository. Therefore, completing the entire upstream pipeline requires users to independently obtain and reconstruct the required source-data inputs.

Because DHS microdata are subject to data-use restrictions, users must obtain access directly from The DHS Program. These files should not be redistributed through this repository.
```