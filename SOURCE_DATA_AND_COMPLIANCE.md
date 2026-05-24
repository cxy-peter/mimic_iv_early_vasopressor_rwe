# Source data and compliance: MIMIC-IV project

MIMIC-IV is credentialed-access data. Do **not** upload raw MIMIC-IV tables or row-level patient datasets to GitHub.

Required local source files:

```text
mimic-iv-3.1/
├── hosp/
│   ├── admissions.csv.gz
│   ├── patients.csv.gz
│   ├── diagnoses_icd.csv.gz
│   └── labevents.csv.gz
└── icu/
    ├── icustays.csv.gz
    ├── inputevents.csv.gz
    ├── d_items.csv.gz
    └── chartevents.csv.gz
```

The project package includes only:

- code,
- notebooks,
- source-data manifest,
- aggregate CSV summaries,
- figures and reports.

The `.gitignore` file blocks raw MIMIC tables and row-level derived cohorts, including `mimic_icu_vasopressor_cohort.csv` and `analysis_dataset_with_weights.csv`.
