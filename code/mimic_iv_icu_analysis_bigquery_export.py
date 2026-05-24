#!/usr/bin/env python3
"""MIMIC-IV ICU treatment strategy analysis.
Run after exporting the BigQuery cohort to CSV.

Example:
python mimic_iv_icu_analysis.py --input mimic_icu_vasopressor_cohort.csv --outdir results/mimic_iv --bootstrap 200
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm


def smd_numeric(x, t, w=None):
    if w is None:
        w = np.ones(len(x))
    mask1 = (t == 1) & np.isfinite(x)
    mask0 = (t == 0) & np.isfinite(x)
    if mask1.sum() < 2 or mask0.sum() < 2:
        return np.nan
    def wm(v, ww): return np.sum(v * ww) / np.sum(ww)
    m1, m0 = wm(x[mask1], w[mask1]), wm(x[mask0], w[mask0])
    v1 = wm((x[mask1]-m1)**2, w[mask1])
    v0 = wm((x[mask0]-m0)**2, w[mask0])
    return (m1 - m0) / np.sqrt((v1 + v0) / 2 + 1e-12)


def make_design(df, covariates):
    cat_cols = [c for c in covariates if df[c].dtype == 'object' or str(df[c].dtype).startswith('category')]
    num_cols = [c for c in covariates if c not in cat_cols]
    pre = ColumnTransformer([
        ('num', Pipeline([('imp', SimpleImputer(strategy='median')), ('scaler', StandardScaler())]), num_cols),
        ('cat', Pipeline([('imp', SimpleImputer(strategy='most_frequent')), ('oh', OneHotEncoder(handle_unknown='ignore'))]), cat_cols),
    ])
    return pre, num_cols, cat_cols


def aipw_binary(df, y_col, t_col, covariates, seed=2026):
    y = df[y_col].astype(float).to_numpy()
    t = df[t_col].astype(int).to_numpy()
    pre, _, _ = make_design(df, covariates)
    ps_model = Pipeline([('pre', pre), ('lr', LogisticRegression(max_iter=2000, C=1.0))])
    ps_model.fit(df[covariates], t)
    e = np.clip(ps_model.predict_proba(df[covariates])[:, 1], 0.01, 0.99)

    outcome_pipe = lambda: Pipeline([('pre', make_design(df, covariates)[0]), ('lr', LogisticRegression(max_iter=2000, C=1.0))])
    m1_model = outcome_pipe(); m0_model = outcome_pipe()
    m1_model.fit(df.loc[t == 1, covariates], y[t == 1])
    m0_model.fit(df.loc[t == 0, covariates], y[t == 0])
    m1 = m1_model.predict_proba(df[covariates])[:, 1]
    m0 = m0_model.predict_proba(df[covariates])[:, 1]
    psi = (m1 - m0) + t * (y - m1) / e - (1 - t) * (y - m0) / (1 - e)
    return psi.mean(), e, m1, m0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', default='results/mimic_iv')
    ap.add_argument('--bootstrap', type=int, default=200)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input)

    required = ['early_vasopressor', 'hospital_expire_flag', 'age', 'sex']
    for col in required:
        if col not in df.columns:
            raise ValueError(f'Missing required column: {col}')

    # Basic cleaning
    df = df.copy()
    df['early_vasopressor'] = df['early_vasopressor'].astype(int)
    df['hospital_expire_flag'] = df['hospital_expire_flag'].astype(int)

    covariates = [c for c in [
        'age','sex','admission_type','charlson','sofa_24h','hr_mean','mbp_mean','rr_mean','temp_mean',
        'creatinine_mean','lactate_mean','wbc_mean','platelet_mean','mechanical_vent','icu_los_days'
    ] if c in df.columns]

    # Cohort flow / descriptive counts
    flow = pd.DataFrame({
        'metric': ['analytic_n','treated_n','control_n','treated_rate','overall_mortality','treated_mortality','control_mortality'],
        'value': [
            len(df), int(df['early_vasopressor'].sum()), int((1-df['early_vasopressor']).sum()),
            df['early_vasopressor'].mean(), df['hospital_expire_flag'].mean(),
            df.loc[df.early_vasopressor==1, 'hospital_expire_flag'].mean(),
            df.loc[df.early_vasopressor==0, 'hospital_expire_flag'].mean()
        ]
    })
    flow.to_csv(outdir/'cohort_flow.csv', index=False)

    pre, _, _ = make_design(df, covariates)
    ps_model = Pipeline([('pre', pre), ('lr', LogisticRegression(max_iter=2000, C=0.5))])
    ps_model.fit(df[covariates], df['early_vasopressor'])
    ps = np.clip(ps_model.predict_proba(df[covariates])[:,1], 0.01, 0.99)
    df['ps'] = ps
    p_treat = df['early_vasopressor'].mean()
    df['iptw'] = np.where(df['early_vasopressor']==1, p_treat/ps, (1-p_treat)/(1-ps))
    lo, hi = df['iptw'].quantile([0.01, 0.99])
    df['iptw_trim'] = df['iptw'].clip(lo, hi)

    # SMD before/after
    smd_rows = []
    for c in covariates:
        if pd.api.types.is_numeric_dtype(df[c]):
            smd_rows.append({'variable': c, 'smd_before': smd_numeric(df[c].to_numpy(float), df.early_vasopressor.to_numpy()),
                             'smd_after': smd_numeric(df[c].to_numpy(float), df.early_vasopressor.to_numpy(), df.iptw_trim.to_numpy())})
    pd.DataFrame(smd_rows).to_csv(outdir/'table1_smd_before_after.csv', index=False)

    # Propensity overlap plot
    plt.figure(figsize=(7,4))
    plt.hist(df.loc[df.early_vasopressor==1, 'ps'], bins=30, alpha=0.5, density=True, label='Early VP')
    plt.hist(df.loc[df.early_vasopressor==0, 'ps'], bins=30, alpha=0.5, density=True, label='No early VP')
    plt.xlabel('Propensity score'); plt.ylabel('Density'); plt.title('Propensity score overlap'); plt.legend(); plt.tight_layout()
    plt.savefig(outdir/'propensity_overlap.png', dpi=200); plt.close()

    # Love plot
    smd_df = pd.DataFrame(smd_rows).dropna()
    if len(smd_df):
        smd_df = smd_df.assign(max_abs=lambda d: np.maximum(d.smd_before.abs(), d.smd_after.abs())).sort_values('max_abs')
        plt.figure(figsize=(7, max(4, 0.25*len(smd_df))))
        y = np.arange(len(smd_df))
        plt.scatter(smd_df.smd_before, y, label='Before')
        plt.scatter(smd_df.smd_after, y, label='After')
        plt.axvline(0.1, linestyle='--'); plt.axvline(-0.1, linestyle='--')
        plt.yticks(y, smd_df.variable); plt.xlabel('Standardized mean difference'); plt.legend(); plt.tight_layout()
        plt.savefig(outdir/'love_plot.png', dpi=200); plt.close()

    # Weighted outcome model
    X = pd.get_dummies(df[['early_vasopressor'] + covariates], drop_first=True)
    X = X.apply(pd.to_numeric, errors='coerce').fillna(X.median(numeric_only=True))
    X = sm.add_constant(X, has_constant='add')
    y = df['hospital_expire_flag']
    glm = sm.GLM(y, X, family=sm.families.Binomial(), var_weights=df['iptw_trim'])
    res = glm.fit(cov_type='HC3')
    pd.DataFrame({'term': res.params.index, 'coef': res.params.values, 'se': res.bse.values,
                  'or': np.exp(res.params.values), 'p': res.pvalues.values}).to_csv(outdir/'weighted_logistic_results.csv', index=False)

    # AIPW risk difference
    ate, e, m1, m0 = aipw_binary(df, 'hospital_expire_flag', 'early_vasopressor', covariates)
    boot = []
    rng = np.random.default_rng(2026)
    for b in range(args.bootstrap):
        idx = rng.integers(0, len(df), len(df))
        try:
            boot.append(aipw_binary(df.iloc[idx].reset_index(drop=True), 'hospital_expire_flag', 'early_vasopressor', covariates)[0])
        except Exception:
            pass
    lo_ci, hi_ci = (np.percentile(boot, [2.5, 97.5]) if len(boot) > 30 else [np.nan, np.nan])
    pd.DataFrame([{'estimand': 'AIPW risk difference: early VP vs no early VP', 'estimate': ate, 'ci_low': lo_ci, 'ci_high': hi_ci,
                   'bootstrap_success_n': len(boot)}]).to_csv(outdir/'aipw_results.csv', index=False)

    print(f'Done. Outputs saved to {outdir}')

if __name__ == '__main__':
    main()
