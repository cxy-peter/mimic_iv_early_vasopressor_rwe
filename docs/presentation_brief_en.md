# Presentation Brief - MIMIC-IV ICU Early Vasopressor RWE

## One-sentence takeaway
This project compares ICU patients who started vasopressors within 6 hours of ICU admission with those who did not start within 6 hours; the crude mortality difference is large, but adjusted analyses suggest it mainly reflects baseline severity rather than a clear treatment harm signal.

## Why this matters
In ICU RWE, treatment timing is strongly tied to clinical severity. Without defining time zero and a clear treatment window, an analysis can easily misread sicker patients' worse outcomes as treatment harm.

## Groups compared

| Group | Definition | N | Crude mortality |
|---|---|---:|---:|
| Early vasopressor | First vasopressor within 0-6h after ICU admission | 3,772 | 32.9% |
| Comparator | No vasopressor initiation within 0-6h after ICU admission | 5,730 | 25.2% |

Important: the comparator is not a never-treated group. It means no initiation before the 6-hour landmark.

## Analysis workflow
1. Build adult first-ICU sepsis-coded cohort.
2. Set ICU admission as time zero.
3. Define early vasopressor initiation in the 0-6h window.
4. Extract baseline clinical features.
5. Estimate propensity scores.
6. Check SMD balance before and after IPTW.
7. Estimate weighted logistic/AIPW risk difference.
8. Interpret crude and adjusted results separately.

## Key results
- Analytic cohort: 9,502 ICU stays.
- Early vasopressor group: 3,772; comparator: 5,730.
- Crude mortality: 32.9% vs 25.2%.
- Balance improved: max absolute SMD 0.288 before IPTW to 0.046 after IPTW.
- AIPW risk difference: 1.17 pp, 95% CI -0.52 to 2.87.

## Management interpretation
The project is not trying to prove whether vasopressors work. It demonstrates how an RWE analyst should structure an ICU treatment-strategy question: define time zero, define the treatment window, create a credible comparator, check balance, and avoid interpreting crude mortality differences as causal.

## Suggested resume line
Built a target-trial-style MIMIC-IV ICU RWE analysis comparing early vasopressor initiation within 0-6h after ICU admission versus no initiation by 6h; improved covariate balance after IPTW (max |SMD| 0.288 to 0.046) and estimated adjusted AIPW RD of 1.17 pp (95% CI -0.52 to 2.87).
