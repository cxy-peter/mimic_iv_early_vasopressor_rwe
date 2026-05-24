# MIMIC-IV source data manifest

## Raw data status

Raw MIMIC-IV files are not included in this repository because they require PhysioNet credentialed access and are governed by a data use agreement.

## Required local folder structure

After credentialed download, the project expects:

```text
/path/to/mimic-iv-3.1/
  hosp/admissions.csv.gz
  hosp/patients.csv.gz
  hosp/diagnoses_icd.csv.gz
  hosp/labevents.csv.gz
  icu/icustays.csv.gz
  icu/inputevents.csv.gz
  icu/d_items.csv.gz
  icu/chartevents.csv.gz
```

## Source tables used

| MIMIC module | Tables | Role |
|---|---|---|
| `hosp` | `admissions`, `patients`, `diagnoses_icd`, `labevents` | demographics, outcome, diagnosis proxy, labs |
| `icu` | `icustays`, `inputevents`, `d_items`, `chartevents` | ICU stays, vasopressor exposure, vitals |

## Derived variables

- `early_vasopressor`: first vasopressor within 0--6 hours after ICU admission.
- `hospital_expire_flag`: in-hospital mortality.
- First-6h vitals/labs: HR, MAP/MBP, RR, temperature, SpO2, lactate, creatinine, WBC, platelet, bilirubin, albumin.

## Compliance note

Only aggregate outputs and scripts should be pushed to GitHub. Patient-level files such as `mimic_icu_vasopressor_cohort.csv` and `analysis_dataset_with_weights.csv` are intentionally ignored by `.gitignore`.
