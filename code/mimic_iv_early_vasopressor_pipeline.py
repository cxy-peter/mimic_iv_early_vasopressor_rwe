#!/usr/bin/env python3
"""
MIMIC-IV local CSV pipeline for a resume-ready RWE / target-trial-emulation project.

Project: Early vasopressor initiation within 0-6 hours of ICU admission and in-hospital mortality
among adult, first-ICU, sepsis-coded MIMIC-IV patients.

Groups compared:
  Treated:    vasopressor started within 0-6 hours after ICU admission.
  Comparator: no vasopressor started within 0-6 hours after ICU admission.
              This does not necessarily mean "never vasopressor"; delayed vasopressor use after 6h
              can be handled in sensitivity analyses.

Important caveat:
  This is a portfolio / resume project template, not a publishable clinical study without clinical validation.
  The local fallback uses diagnosis-code sepsis proxy and ICU admission as time zero. For a stronger version,
  use mimic-code derived concepts such as sepsis3, first_day_sofa, charlson, ventilation, vitalsign, chemistry.

Expected folder structure:
  /path/to/mimiciv/
    hosp/admissions.csv.gz
    hosp/patients.csv.gz
    hosp/diagnoses_icd.csv.gz
    hosp/labevents.csv.gz
    icu/icustays.csv.gz
    icu/inputevents.csv.gz
    icu/d_items.csv.gz
    icu/chartevents.csv.gz

Run:
  python code/mimic_iv_early_vasopressor_pipeline.py --mimic_root /path/to/mimic-iv-3.1 --outdir outputs_reproduced --treatment_window_hours 6

Analysis-only if cohort CSV already exists:
  python code/mimic_iv_early_vasopressor_pipeline.py --cohort_csv outputs_reproduced/mimic_icu_vasopressor_cohort.csv --outdir outputs_reproduced --analysis_only
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def sql_path(path: str | Path) -> str:
    return str(path).replace('\\', '/').replace("'", "''")


def resolve_mimic_root(root: Path) -> Path:
    root = Path(root)
    if (root / 'hosp').exists() and (root / 'icu').exists():
        return root
    if root.exists():
        for p in root.rglob('hosp'):
            candidate = p.parent
            if (candidate / 'icu').exists():
                return candidate
    return root


def find_table(root: Path, module: str, table: str) -> str:
    for suffix in ['.csv.gz', '.csv']:
        p = root / module / f'{table}{suffix}'
        if p.exists():
            return str(p)
    raise FileNotFoundError(f'Missing {module}/{table}.csv.gz or .csv under {root}')


def check_required_files(root: Path) -> dict[str, str]:
    required = {
        'icustays': ('icu', 'icustays'),
        'inputevents': ('icu', 'inputevents'),
        'd_items': ('icu', 'd_items'),
        'chartevents': ('icu', 'chartevents'),
        'admissions': ('hosp', 'admissions'),
        'patients': ('hosp', 'patients'),
        'diagnoses_icd': ('hosp', 'diagnoses_icd'),
        'labevents': ('hosp', 'labevents'),
    }
    return {name: find_table(root, module, table) for name, (module, table) in required.items()}


def create_views(con: duckdb.DuckDBPyConnection, paths: dict[str, str]) -> None:
    """Create DuckDB views over source CSV/CSV.GZ files.

    v2 uses a tolerant reader because Windows-extracted large CSV files can contain a small number
    of malformed physical lines. For a publishable study, re-download the official files and audit
    skipped rows; for this portfolio pipeline, ignore_errors/null_padding prevents one broken line
    from stopping the full run.
    """
    csv_options = 'union_by_name=true, ignore_errors=true, null_padding=true'
    for view_name, file_path in paths.items():
        print(f'Creating view {view_name}: {file_path}')
        con.execute(
            f"CREATE OR REPLACE VIEW {view_name} AS "
            f"SELECT * FROM read_csv_auto('{sql_path(file_path)}', {csv_options});"
        )


def build_source_cohort(con: duckdb.DuckDBPyConnection, outdir: Path, treatment_window_hours: int) -> pd.DataFrame:
    """Build adult first-ICU sepsis-coded cohort and derive early vasopressor exposure."""
    con.execute(f"""
    CREATE OR REPLACE TEMP TABLE first_icu AS
    SELECT * EXCLUDE(rn)
    FROM (
      SELECT i.*, ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY intime) AS rn
      FROM icustays i
    )
    WHERE rn = 1;

    CREATE OR REPLACE TEMP TABLE sepsis_hadm AS
    SELECT DISTINCT hadm_id
    FROM diagnoses_icd
    WHERE
      (icd_version = 9 AND (icd_code LIKE '038%' OR icd_code IN ('99591','99592','78552')))
      OR
      (icd_version = 10 AND (icd_code LIKE 'A40%' OR icd_code LIKE 'A41%' OR icd_code LIKE 'R652%'));

    CREATE OR REPLACE TEMP TABLE base_cohort AS
    SELECT
      i.subject_id,
      i.hadm_id,
      i.stay_id,
      i.first_careunit,
      i.intime,
      i.outtime,
      a.admittime,
      a.dischtime,
      a.deathtime,
      a.hospital_expire_flag,
      a.admission_type,
      p.gender AS sex,
      p.anchor_age AS age,
      date_diff('hour', i.intime, i.outtime) AS icu_los_hours,
      date_diff('hour', i.intime, a.dischtime) AS hosp_los_from_icu_hours,
      date_diff('hour', a.admittime, a.dischtime) / 24.0 AS hosp_los_days,
      COUNT(d.icd_code) AS diagnosis_count
    FROM first_icu i
    JOIN admissions a USING (subject_id, hadm_id)
    JOIN patients p USING (subject_id)
    JOIN sepsis_hadm s USING (hadm_id)
    LEFT JOIN diagnoses_icd d USING (subject_id, hadm_id)
    WHERE p.anchor_age >= 18
    GROUP BY
      i.subject_id, i.hadm_id, i.stay_id, i.first_careunit, i.intime, i.outtime,
      a.admittime, a.dischtime, a.deathtime, a.hospital_expire_flag, a.admission_type,
      p.gender, p.anchor_age;
    """)

    con.execute("""
    CREATE OR REPLACE TEMP TABLE vaso_items AS
    SELECT DISTINCT itemid, label
    FROM d_items
    WHERE lower(label) LIKE '%norepinephrine%'
       OR lower(label) LIKE '%epinephrine%'
       OR lower(label) LIKE '%vasopressin%'
       OR lower(label) LIKE '%phenylephrine%'
       OR lower(label) LIKE '%dopamine%';

    CREATE OR REPLACE TEMP TABLE first_vaso AS
    SELECT
      b.stay_id,
      MIN(ie.starttime) AS first_vp_time,
      MIN(date_diff('minute', b.intime, ie.starttime)) / 60.0 AS first_vp_hours
    FROM base_cohort b
    JOIN inputevents ie ON b.stay_id = ie.stay_id
    JOIN vaso_items vi ON ie.itemid = vi.itemid
    WHERE ie.starttime >= b.intime
      AND ie.starttime < b.intime + INTERVAL 24 HOUR
      AND COALESCE(ie.amount, 0) > 0
    GROUP BY b.stay_id;
    """)

    con.execute(f"""
    CREATE OR REPLACE TEMP TABLE vital_6h AS
    SELECT
      b.stay_id,
      AVG(CASE WHEN ce.itemid = 220045 THEN ce.valuenum END) AS hr_mean,
      AVG(CASE WHEN ce.itemid IN (220052, 220181) THEN ce.valuenum END) AS mbp_mean,
      AVG(CASE WHEN ce.itemid = 220210 THEN ce.valuenum END) AS rr_mean,
      AVG(CASE WHEN ce.itemid = 220277 THEN ce.valuenum END) AS spo2_mean,
      AVG(CASE WHEN ce.itemid = 223762 THEN ce.valuenum END) AS temp_c_mean,
      AVG(CASE WHEN ce.itemid = 223761 THEN (ce.valuenum - 32) * 5.0 / 9.0 END) AS temp_f_as_c_mean
    FROM base_cohort b
    LEFT JOIN chartevents ce
      ON b.stay_id = ce.stay_id
     AND ce.charttime >= b.intime
     AND ce.charttime < b.intime + INTERVAL {treatment_window_hours} HOUR
     AND ce.valuenum IS NOT NULL
    GROUP BY b.stay_id;

    CREATE OR REPLACE TEMP TABLE lab_6h AS
    SELECT
      b.stay_id,
      AVG(CASE WHEN le.itemid = 50813 THEN le.valuenum END) AS lactate_mean,
      AVG(CASE WHEN le.itemid = 50912 THEN le.valuenum END) AS creatinine_mean,
      AVG(CASE WHEN le.itemid IN (51300, 51301) THEN le.valuenum END) AS wbc_mean,
      AVG(CASE WHEN le.itemid = 51265 THEN le.valuenum END) AS platelet_mean,
      AVG(CASE WHEN le.itemid = 50885 THEN le.valuenum END) AS bilirubin_total_mean,
      AVG(CASE WHEN le.itemid = 50862 THEN le.valuenum END) AS albumin_mean
    FROM base_cohort b
    LEFT JOIN labevents le
      ON b.hadm_id = le.hadm_id
     AND le.charttime >= b.intime
     AND le.charttime < b.intime + INTERVAL {treatment_window_hours} HOUR
     AND le.valuenum IS NOT NULL
    GROUP BY b.stay_id;
    """)

    df = con.execute(f"""
    SELECT
      b.subject_id,
      b.hadm_id,
      b.stay_id,
      b.age,
      b.sex,
      b.admission_type,
      b.first_careunit,
      b.intime,
      b.outtime,
      b.admittime,
      b.dischtime,
      b.deathtime,
      b.hospital_expire_flag,
      b.icu_los_hours / 24.0 AS icu_los_days,
      b.hosp_los_days,
      b.diagnosis_count,
      fv.first_vp_time,
      fv.first_vp_hours,
      CASE WHEN fv.first_vp_hours >= 0 AND fv.first_vp_hours < {treatment_window_hours} THEN 1 ELSE 0 END AS early_vasopressor,
      CASE WHEN fv.first_vp_hours >= {treatment_window_hours} AND fv.first_vp_hours < 24 THEN 1 ELSE 0 END AS delayed_vasopressor_6_24h,
      CASE WHEN fv.first_vp_hours >= 0 AND fv.first_vp_hours < 24 THEN 1 ELSE 0 END AS received_vasopressor_24h,
      v.hr_mean,
      v.mbp_mean,
      v.rr_mean,
      COALESCE(v.temp_c_mean, v.temp_f_as_c_mean) AS temp_mean,
      v.spo2_mean,
      l.creatinine_mean,
      l.lactate_mean,
      l.wbc_mean,
      l.platelet_mean,
      l.bilirubin_total_mean,
      l.albumin_mean
    FROM base_cohort b
    LEFT JOIN first_vaso fv USING (stay_id)
    LEFT JOIN vital_6h v USING (stay_id)
    LEFT JOIN lab_6h l USING (stay_id)
    WHERE b.icu_los_hours >= {treatment_window_hours}
      AND b.hosp_los_from_icu_hours >= {treatment_window_hours};
    """).fetchdf()

    df['early_vasopressor'] = df['early_vasopressor'].astype(int)
    df['hospital_expire_flag'] = df['hospital_expire_flag'].astype(int)
    df.to_csv(outdir / 'mimic_icu_vasopressor_cohort.csv', index=False)

    con.execute('SELECT * FROM vaso_items ORDER BY label').fetchdf().to_csv(outdir / 'vasopressor_itemids_detected.csv', index=False)
    con.execute('SELECT * FROM base_cohort').fetchdf().to_csv(outdir / '01_base_sepsis_first_icu_cohort.csv', index=False)
    con.execute('SELECT * FROM first_vaso').fetchdf().to_csv(outdir / '02_first_vasopressor_by_stay.csv', index=False)
    con.execute('SELECT * FROM vital_6h').fetchdf().to_csv(outdir / '03_vitals_6h.csv', index=False)
    con.execute('SELECT * FROM lab_6h').fetchdf().to_csv(outdir / '03_labs_6h.csv', index=False)
    return df


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(w)
    if ok.sum() == 0 or w[ok].sum() == 0:
        return np.nan
    return float(np.sum(x[ok] * w[ok]) / np.sum(w[ok]))


def smd_numeric(x: np.ndarray, t: np.ndarray, w: np.ndarray | None = None) -> float:
    if w is None:
        w = np.ones(len(x))
    x = np.asarray(x, dtype=float)
    t = np.asarray(t, dtype=int)
    w = np.asarray(w, dtype=float)
    m1 = (t == 1) & np.isfinite(x)
    m0 = (t == 0) & np.isfinite(x)
    if m1.sum() < 2 or m0.sum() < 2:
        return np.nan
    mu1, mu0 = weighted_mean(x[m1], w[m1]), weighted_mean(x[m0], w[m0])
    v1 = weighted_mean((x[m1] - mu1) ** 2, w[m1])
    v0 = weighted_mean((x[m0] - mu0) ** 2, w[m0])
    return float((mu1 - mu0) / math.sqrt((v1 + v0) / 2.0 + 1e-12))


def get_covariates(df: pd.DataFrame) -> list[str]:
    # Baseline/pre-exposure-style variables only. Avoid LOS as adjustment variable.
    candidates = [
        'age', 'sex', 'admission_type', 'first_careunit', 'diagnosis_count',
        'hr_mean', 'mbp_mean', 'rr_mean', 'temp_mean', 'spo2_mean',
        'creatinine_mean', 'lactate_mean', 'wbc_mean', 'platelet_mean',
        'bilirubin_total_mean', 'albumin_mean',
    ]
    return [c for c in candidates if c in df.columns]


def make_preprocessor(df: pd.DataFrame, covariates: Iterable[str]) -> ColumnTransformer:
    covariates = list(covariates)
    cat_cols = [c for c in covariates if df[c].dtype == 'object' or str(df[c].dtype).startswith('category')]
    num_cols = [c for c in covariates if c not in cat_cols]
    return ColumnTransformer([
        ('num', Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler())]), num_cols),
        ('cat', Pipeline([('imp', SimpleImputer(strategy='most_frequent')), ('oh', OneHotEncoder(handle_unknown='ignore'))]), cat_cols),
    ])


def run_analysis(df: pd.DataFrame, outdir: Path, treatment_col: str = 'early_vasopressor', outcome_col: str = 'hospital_expire_flag') -> None:
    df = df.copy()
    covariates = get_covariates(df)
    required = [treatment_col, outcome_col, 'age', 'sex']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f'Missing required columns: {missing}')
    if df[treatment_col].nunique() < 2:
        raise ValueError('Treatment has only one level. Check cohort definition or treatment window.')
    if not covariates:
        raise ValueError('No baseline covariates found.')

    flow = pd.DataFrame({
        'metric': ['analytic_n', 'treated_n', 'control_n', 'treated_rate', 'overall_mortality', 'treated_mortality', 'control_mortality'],
        'value': [len(df), int(df[treatment_col].sum()), int((1 - df[treatment_col]).sum()), float(df[treatment_col].mean()), float(df[outcome_col].mean()), float(df.loc[df[treatment_col] == 1, outcome_col].mean()), float(df.loc[df[treatment_col] == 0, outcome_col].mean())],
    })
    flow.to_csv(outdir / 'cohort_flow.csv', index=False)
    flow.to_csv(outdir / 'cohort_flow_summary.csv', index=False)
    crude = df.groupby(treatment_col)[outcome_col].agg(['count', 'mean']).reset_index()
    crude.to_csv(outdir / 'crude_mortality_by_treatment.csv', index=False)
    crude_summary = crude.rename(columns={treatment_col: 'early_vasopressor_0_6h', 'mean': 'hospital_mortality_rate'})
    crude_summary['group_label'] = crude_summary['early_vasopressor_0_6h'].map({0: 'No early vasopressor initiation within 0-6h', 1: 'Early vasopressor initiation within 0-6h'})
    crude_summary = crude_summary[['early_vasopressor_0_6h', 'group_label', 'count', 'hospital_mortality_rate']].rename(columns={'count': 'n'})
    crude_summary.to_csv(outdir / 'crude_mortality_by_treatment_summary.csv', index=False)

    pre = make_preprocessor(df, covariates)
    ps_model = Pipeline([('pre', pre), ('lr', LogisticRegression(max_iter=3000, C=0.5, class_weight='balanced'))])
    ps_model.fit(df[covariates], df[treatment_col].astype(int))
    ps = np.clip(ps_model.predict_proba(df[covariates])[:, 1], 0.01, 0.99)
    df['propensity_score'] = ps
    p_treat = df[treatment_col].mean()
    df['iptw_stabilized'] = np.where(df[treatment_col] == 1, p_treat / ps, (1 - p_treat) / (1 - ps))
    lo, hi = df['iptw_stabilized'].quantile([0.01, 0.99])
    df['iptw_trimmed'] = df['iptw_stabilized'].clip(lo, hi)
    df.to_csv(outdir / 'analysis_dataset_with_weights.csv', index=False)

    smd_rows = []
    for c in covariates:
        if pd.api.types.is_numeric_dtype(df[c]):
            smd_rows.append({'variable': c, 'smd_before': smd_numeric(df[c].to_numpy(float), df[treatment_col].to_numpy(int)), 'smd_after_iptw_trimmed': smd_numeric(df[c].to_numpy(float), df[treatment_col].to_numpy(int), df['iptw_trimmed'].to_numpy(float))})
    smd_df = pd.DataFrame(smd_rows).sort_values('variable')
    smd_df.to_csv(outdir / 'table1_smd_before_after.csv', index=False)
    if not smd_df.empty:
        balance_summary = pd.DataFrame({
            'balance_metric': ['numeric covariates with |SMD| < 0.1 after IPTW', 'max |SMD| before IPTW', 'max |SMD| after IPTW'],
            'value': [f"{int((smd_df['smd_after_iptw_trimmed'].abs() < 0.1).sum())}/{int(smd_df['smd_after_iptw_trimmed'].notna().sum())}", round(float(smd_df['smd_before'].abs().max()), 3), round(float(smd_df['smd_after_iptw_trimmed'].abs().max()), 3)]
        })
        balance_summary.to_csv(outdir / 'balance_summary.csv', index=False)

    plt.figure(figsize=(7, 4))
    plt.hist(df.loc[df[treatment_col] == 1, 'propensity_score'], bins=30, alpha=0.5, density=True, label='Early vasopressor')
    plt.hist(df.loc[df[treatment_col] == 0, 'propensity_score'], bins=30, alpha=0.5, density=True, label='No early vasopressor')
    plt.xlabel('Propensity score')
    plt.ylabel('Density')
    plt.title('Propensity score overlap')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / 'propensity_overlap.png', dpi=200)
    plt.close()

    if not smd_df.empty:
        plot_df = smd_df.dropna().assign(max_abs=lambda d: np.maximum(d['smd_before'].abs(), d['smd_after_iptw_trimmed'].abs())).sort_values('max_abs')
        plt.figure(figsize=(7, max(4, 0.30 * len(plot_df))))
        y = np.arange(len(plot_df))
        plt.scatter(plot_df['smd_before'], y, label='Before')
        plt.scatter(plot_df['smd_after_iptw_trimmed'], y, label='After IPTW')
        plt.axvline(0.1, linestyle='--')
        plt.axvline(-0.1, linestyle='--')
        plt.yticks(y, plot_df['variable'])
        plt.xlabel('Standardized mean difference')
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / 'love_plot.png', dpi=200)
        plt.close()

    # Weighted logistic regression. v2/v3 fix: convert dummy-coded design matrix to float.
    weighted_logistic_path = outdir / 'weighted_logistic_results.csv'
    try:
        X = pd.get_dummies(df[[treatment_col] + covariates], drop_first=True)
        X = X.apply(pd.to_numeric, errors='coerce')
        X = X.fillna(X.median(numeric_only=True)).astype(float)
        X = sm.add_constant(X, has_constant='add')
        y = df[outcome_col].astype(int)
        glm = sm.GLM(y, X, family=sm.families.Binomial(), var_weights=df['iptw_trimmed'].astype(float))
        res = glm.fit(cov_type='HC3')
        weighted_logistic = pd.DataFrame({'term': res.params.index, 'coef': res.params.values, 'se': res.bse.values, 'or': np.exp(res.params.values), 'p': res.pvalues.values})
        weighted_logistic.to_csv(weighted_logistic_path, index=False)
    except Exception as e:
        pd.DataFrame([{'error': repr(e), 'note': 'Weighted logistic model failed, but AIPW will continue.'}]).to_csv(weighted_logistic_path, index=False)
        print(f'WARNING: weighted logistic failed but AIPW will continue: {e}')

    A = df[treatment_col].astype(int).to_numpy()
    Y = df[outcome_col].astype(int).to_numpy()
    X_cov = df[covariates].copy()
    outcome_data = X_cov.copy()
    outcome_data[treatment_col] = A
    outcome_covariates = covariates + [treatment_col]
    q_pre = make_preprocessor(outcome_data, outcome_covariates)
    q_model = Pipeline([('pre', q_pre), ('gb', GradientBoostingClassifier(random_state=2026))])
    q_model.fit(outcome_data[outcome_covariates], Y)
    x1 = X_cov.copy(); x1[treatment_col] = 1
    x0 = X_cov.copy(); x0[treatment_col] = 0
    m1 = q_model.predict_proba(x1[outcome_covariates])[:, 1]
    m0 = q_model.predict_proba(x0[outcome_covariates])[:, 1]
    pseudo = m1 - m0 + A / ps * (Y - m1) - (1 - A) / (1 - ps) * (Y - m0)
    rd = float(np.mean(pseudo))
    se = float(np.std(pseudo, ddof=1) / math.sqrt(len(df)))

    aipw_df = pd.DataFrame([{'estimand': 'AIPW risk difference: early vasopressor vs no early vasopressor', 'estimate': rd, 'ci_low': rd - 1.96 * se, 'ci_high': rd + 1.96 * se, 'n': len(df), 'treated_n': int(A.sum()), 'control_n': int((1 - A).sum())}])
    aipw_df.to_csv(outdir / 'aipw_results.csv', index=False)
    aipw_summary = aipw_df.copy()
    for c in ['estimate', 'ci_low', 'ci_high']:
        aipw_summary[c + '_pp'] = aipw_summary[c] * 100
    aipw_summary = aipw_summary[['estimand', 'estimate_pp', 'ci_low_pp', 'ci_high_pp', 'n', 'treated_n', 'control_n']]
    aipw_summary.to_csv(outdir / 'aipw_results_summary.csv', index=False)

    treated_mort = float(df.loc[df[treatment_col] == 1, outcome_col].mean())
    control_mort = float(df.loc[df[treatment_col] == 0, outcome_col].mean())
    with open(outdir / 'resume_numbers.md', 'w', encoding='utf-8') as f:
        f.write('# Resume-ready numbers\n\n')
        f.write(f'- Analytic cohort N = {len(df):,}\n')
        f.write(f'- Early vasopressor group N = {int(A.sum()):,}\n')
        f.write(f'- No early vasopressor group N = {int((1 - A).sum()):,}\n')
        f.write(f'- Crude hospital mortality: early = {treated_mort:.3f}, control = {control_mort:.3f}\n')
        f.write(f'- AIPW risk difference = {rd:.4f} ({rd - 1.96 * se:.4f}, {rd + 1.96 * se:.4f})\n')
        f.write('\nSuggested wording: observational association, not causal proof; confounding by indication remains possible.\n')



def main() -> None:
    parser = argparse.ArgumentParser(description='Unified extraction + analysis pipeline for the MIMIC-IV early vasopressor RWE project.')
    parser.add_argument('--mimic_root', default='', help='Root folder containing MIMIC-IV hosp/ and icu/ folders. Required unless --analysis_only with --cohort_csv is used.')
    parser.add_argument('--cohort_csv', default='', help='Optional existing mimic_icu_vasopressor_cohort.csv. Use with --analysis_only to rerun downstream analysis only.')
    parser.add_argument('--outdir', default='outputs_reproduced', help='Output folder for reproduced cohort, diagnostics, and summaries.')
    parser.add_argument('--treatment_window_hours', type=int, default=6, help='Early-treatment landmark window after ICU admission.')
    parser.add_argument('--memory_limit', default='8GB', help='DuckDB memory limit for local extraction.')
    parser.add_argument('--analysis_only', action='store_true', help='Skip raw-data extraction and run only downstream analysis from --cohort_csv.')
    parser.add_argument('--extract_only', action='store_true', help='Extract the cohort but skip propensity-score / AIPW analysis.')
    parser.add_argument('--check_only', action='store_true', help='Only check that required source files exist; do not extract or analyze.')
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)

    if args.analysis_only or args.cohort_csv:
        cohort_path = Path(args.cohort_csv) if args.cohort_csv else (outdir / 'mimic_icu_vasopressor_cohort.csv')
        if not cohort_path.exists():
            raise FileNotFoundError(f'Cannot find cohort CSV for analysis-only mode: {cohort_path}')
        print(f'Analysis-only mode. Reading cohort: {cohort_path}')
        df = pd.read_csv(cohort_path)
        run_analysis(df, outdir)
        print(f'Done. Outputs saved to: {outdir.resolve()}')
        return

    if not args.mimic_root:
        raise ValueError('Provide --mimic_root for extraction mode, or use --analysis_only --cohort_csv for downstream analysis only.')

    root = resolve_mimic_root(Path(args.mimic_root))
    print(f'Using MIMIC root: {root}')
    print(f'Saving outputs to: {outdir}')
    paths = check_required_files(root)
    print('Required files found:')
    for name, path in paths.items():
        print(f'  - {name}: {path}')
    if args.check_only:
        print('Check completed. No extraction or analysis was run.')
        return

    con = duckdb.connect(database=':memory:')
    temp_dir = ensure_dir(outdir / '_duckdb_tmp')
    con.execute(f"PRAGMA temp_directory='{sql_path(str(temp_dir))}';")
    con.execute(f"PRAGMA memory_limit='{args.memory_limit}';")
    create_views(con, paths)
    df = build_source_cohort(con, outdir, args.treatment_window_hours)
    if not args.extract_only:
        run_analysis(df, outdir)
    print(f'Done. Outputs saved to: {outdir.resolve()}')
    print(f'Open this first: {outdir / "resume_numbers.md"}')


if __name__ == '__main__':
    main()
