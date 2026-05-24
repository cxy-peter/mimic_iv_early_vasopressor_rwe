# Run order: MIMIC-IV early vasopressor RWE project

## 0. Install dependencies

```bash
pip install -r environment/requirements_python.txt
```

For the R Markdown notebook:

```r
source("environment/install_r_packages.R")
```

## 1. Check raw MIMIC-IV files

```bash
python code/mimic_iv_early_vasopressor_pipeline.py \
  --mimic_root /path/to/mimic-iv-3.1 \
  --outdir outputs_reproduced \
  --check_only
```

## 2. Full extraction + analysis

```bash
python code/mimic_iv_early_vasopressor_pipeline.py \
  --mimic_root /path/to/mimic-iv-3.1 \
  --outdir outputs_reproduced \
  --treatment_window_hours 6
```

## 3. Analysis-only rerun after cohort extraction

```bash
python code/mimic_iv_early_vasopressor_pipeline.py \
  --cohort_csv outputs_reproduced/mimic_icu_vasopressor_cohort.csv \
  --outdir outputs_reproduced \
  --analysis_only
```

## 4. Render the notebook

Open and knit:

```text
notebooks/01_mimic_iv_early_vasopressor_pipeline.Rmd
```

The notebook documents the target-trial-style design, cohort creation SQL, why the previous split scripts were consolidated, and the final analysis outputs.
