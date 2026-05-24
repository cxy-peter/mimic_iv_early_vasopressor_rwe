# Academic Report - MIMIC-IV ICU Early Vasopressor Treatment Strategy RWE

## Title
Early Vasopressor Initiation and In-hospital Mortality Among Sepsis-coded ICU Stays: A Target-trial-style Real-world Evidence Analysis Using MIMIC-IV

## Abstract
This portfolio project evaluates whether early vasopressor initiation is associated with in-hospital mortality among adult, first-ICU, sepsis-coded stays in MIMIC-IV. ICU admission was treated as time zero. The treatment strategy group was defined as first vasopressor initiation within 0-6 hours after ICU admission. The comparator group was defined as no vasopressor initiation within the same 0-6 hour window. The comparator should not be interpreted as a never-treated group; delayed treatment after the 6-hour landmark can occur. The final analytic cohort contained 9,502 ICU stays, including 3,772 early vasopressor initiators and 5,730 comparators. Crude in-hospital mortality was 32.9% in the early vasopressor group and 25.2% in the comparator group. After propensity-score weighting, all 13 numeric covariates had absolute standardized mean differences below 0.1, with the maximum absolute SMD reduced from 0.288 to 0.046. The AIPW risk difference for early vasopressor initiation versus no initiation by 6 hours was 1.17 percentage points (95% CI: -0.52 to 2.87). The main conclusion is not that vasopressors are harmful; rather, early vasopressor initiation is strongly related to acute illness severity, and crude outcome comparisons should not be interpreted causally without careful design and adjustment.

## Research Question
Among adult first-ICU sepsis-coded stays, what is the adjusted association between early vasopressor initiation within 6 hours of ICU admission and in-hospital mortality?

## Treatment Strategy Contrast

| Component | Definition |
|---|---|
| Time zero | ICU admission time |
| Treatment group | First vasopressor initiation within 0-6 hours after ICU admission |
| Comparator group | No vasopressor initiation within 0-6 hours after ICU admission |
| Outcome | In-hospital mortality |
| Important caveat | Comparator means not treated before the 6-hour landmark, not necessarily never treated |

## Data Source and Cohort
The project uses credentialed MIMIC-IV ICU data. Raw MIMIC-IV tables are not included in this repository. The included code expects a local credentialed MIMIC-IV directory or a BigQuery-derived cohort export. The cohort construction logic restricts to adult first-ICU stays with sepsis coding and sufficient baseline information.

## Methods
The analysis follows a target-trial-style observational workflow:

1. Define eligibility, time zero, treatment strategy, comparator strategy, and outcome.
2. Identify first vasopressor initiation from ICU input event records.
3. Construct baseline features from information available at or before the early ICU window.
4. Estimate propensity scores for early vasopressor initiation.
5. Check covariate balance using standardized mean differences before and after IPTW.
6. Fit propensity-score-weighted logistic regression and estimate an AIPW risk difference.
7. Report crude rates separately from adjusted estimates.

The final unified Python pipeline consolidates earlier extraction and analysis-only scripts. The main fix from the later analysis-only version was retained: after dummy encoding categorical variables, the weighted logistic regression design matrix is explicitly converted to numeric float types before fitting in `statsmodels`. This avoids failures caused by boolean or object dtype columns.

## Results

| Quantity | Result |
|---|---:|
| Analytic cohort | 9,502 ICU stays |
| Early vasopressor group | 3,772 |
| Comparator group | 5,730 |
| Early vasopressor rate | 39.7% |
| Overall in-hospital mortality | 28.3% |
| Crude mortality, early vasopressor | 32.9% |
| Crude mortality, comparator | 25.2% |
| Maximum absolute SMD before IPTW | 0.288 |
| Maximum absolute SMD after IPTW | 0.046 |
| Numeric covariates balanced after IPTW | 13/13 with |SMD| < 0.1 |
| AIPW risk difference | 1.17 pp |
| 95% CI | -0.52 to 2.87 pp |

## Interpretation
Patients requiring early vasopressor therapy are typically more severely ill at baseline. The crude mortality gap should therefore be interpreted as a mixture of treatment timing, severity, and residual differences in clinical status. After design-based adjustment, the estimated risk difference was small and statistically uncertain. The project demonstrates how an RWE analyst should separate crude clinical patterns from adjusted treatment-strategy estimates.

## Limitations
- Sepsis coding and vasopressor exposure definitions are simplified for a portfolio project.
- Residual confounding is likely because treatment initiation is driven by unmeasured clinical severity and bedside judgment.
- Immortal-time and landmark design choices require careful sensitivity analyses.
- The comparator is not a never-treated group.
- Results should not be used as clinical evidence for vasopressor efficacy or safety.

## Reproducibility
See `RUN_ORDER.md` and `notebooks/01_mimic_iv_early_vasopressor_pipeline.Rmd`. Raw MIMIC-IV data are excluded by design; only aggregate outputs and code are included.
