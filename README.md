# MIMIC-IV ICU Treatment Strategy RWE: Early Vasopressor vs No Early Vasopressor by 6 Hours

This project is a portfolio-style real-world evidence (RWE) analysis using MIMIC-IV ICU data. It compares two **treatment strategy groups** among adult, first-ICU, sepsis-coded stays:

| Group | Definition |
|---|---|
| Treatment | First vasopressor initiation within 0--6 hours after ICU admission |
| Comparator | No vasopressor initiation within 0--6 hours after ICU admission |

The comparator is **not** "never received vasopressor." It means the patient did not initiate vasopressor therapy before the 6-hour landmark; delayed initiation after 6 hours can still occur and should be handled as a sensitivity question.

## Why the final code is one script

Earlier drafts had two scripts:

- extraction + analysis pipeline, and
- an analysis-only script.

That was confusing because people could not tell whether v2 or v3 should be run first. The current version consolidates both into one entry point:

```bash
python code/mimic_iv_early_vasopressor_pipeline.py
```

The final script integrates the important fix from the later analysis-only version: the statsmodels weighted logistic design matrix is explicitly converted to numeric/float after dummy encoding. This avoids failures caused by bool/object dtypes in categorical dummy columns.

## Main workflow

1. Check credentialed MIMIC-IV source files.
2. Build the adult first-ICU sepsis-coded cohort.
3. Define ICU admission as time zero.
4. Identify first vasopressor start time from ICU inputevents.
5. Define treatment as initiation within 0--6h.
6. Extract first-6h baseline vitals/labs.
7. Estimate propensity scores.
8. Check balance using SMD before/after IPTW.
9. Fit PS-weighted logistic regression and AIPW risk difference.
10. Export cohort flow, crude rates, balance diagnostics, plots, and resume-ready numbers.

## Run

See [`RUN_ORDER.md`](RUN_ORDER.md). The main notebook is:

```text
notebooks/01_mimic_iv_early_vasopressor_pipeline.Rmd
```

## Outputs included

The `outputs/` folder contains aggregate outputs and generated documents from the completed run. Patient-level MIMIC-IV source data are intentionally not included.

## Resume wording

Use conservative wording: this is an observational treatment-strategy analysis. Early vasopressor use is also a marker of shock severity; crude mortality differences should not be interpreted as vasopressor harm.
