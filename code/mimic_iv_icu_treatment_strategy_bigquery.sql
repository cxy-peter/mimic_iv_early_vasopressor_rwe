-- MIMIC-IV ICU Treatment Strategy Project
-- Early vasopressor initiation within 6h after sepsis onset and in-hospital mortality
-- Replace project/dataset names below before running.
-- Assumes MIT-LCP mimic-code derived tables are available. If not, first build sepsis3, charlson,
-- first_day_sofa, vitalsign, chemistry, vasopressor using mimic-code concept SQL.

DECLARE treatment_window_hours INT64 DEFAULT 6;

-- TODO: replace these dataset names for your environment.
-- Typical BigQuery datasets after access may include:
-- `physionet-data.mimiciv_v3_1_hosp.*`
-- `physionet-data.mimiciv_v3_1_icu.*`
-- Derived tables may be in your own project, e.g. `my_project.mimiciv_derived.*`

WITH base AS (
  SELECT
    icu.subject_id,
    icu.hadm_id,
    icu.stay_id,
    icu.intime,
    icu.outtime,
    DATETIME_DIFF(icu.outtime, icu.intime, HOUR) / 24.0 AS icu_los_days,
    adm.admittime,
    adm.dischtime,
    adm.hospital_expire_flag,
    adm.admission_type,
    pat.anchor_age AS age,
    pat.sex,
    DATETIME_DIFF(adm.dischtime, adm.admittime, HOUR) / 24.0 AS hosp_los_days,
    ROW_NUMBER() OVER (PARTITION BY icu.subject_id ORDER BY icu.intime) AS icu_stay_num
  FROM `physionet-data.mimiciv_v3_1_icu.icustays` icu
  JOIN `physionet-data.mimiciv_v3_1_hosp.admissions` adm
    ON icu.hadm_id = adm.hadm_id
  JOIN `physionet-data.mimiciv_v3_1_hosp.patients` pat
    ON icu.subject_id = pat.subject_id
),
sepsis AS (
  SELECT stay_id, MIN(suspected_infection_time) AS sepsis_time
  FROM `YOUR_PROJECT.mimiciv_derived.sepsis3`
  WHERE sepsis3 = TRUE
  GROUP BY stay_id
),
vp AS (
  SELECT
    stay_id,
    MIN(starttime) AS first_vp_time
  FROM `YOUR_PROJECT.mimiciv_derived.vasopressor`
  GROUP BY stay_id
),
sofa AS (
  SELECT stay_id, sofa AS sofa_24h
  FROM `YOUR_PROJECT.mimiciv_derived.first_day_sofa`
),
charlson AS (
  SELECT subject_id, hadm_id, charlson_comorbidity_index AS charlson
  FROM `YOUR_PROJECT.mimiciv_derived.charlson`
),
vitals AS (
  SELECT
    stay_id,
    AVG(heart_rate) AS hr_mean,
    AVG(mbp) AS mbp_mean,
    AVG(resp_rate) AS rr_mean,
    AVG(temperature) AS temp_mean
  FROM `YOUR_PROJECT.mimiciv_derived.vitalsign`
  GROUP BY stay_id
),
chem AS (
  SELECT
    stay_id,
    AVG(creatinine) AS creatinine_mean,
    AVG(lactate) AS lactate_mean,
    AVG(wbc) AS wbc_mean,
    AVG(platelet) AS platelet_mean
  FROM `YOUR_PROJECT.mimiciv_derived.chemistry`
  GROUP BY stay_id
),
vent AS (
  SELECT
    stay_id,
    MAX(CASE WHEN ventilation_status IS NOT NULL THEN 1 ELSE 0 END) AS mechanical_vent
  FROM `YOUR_PROJECT.mimiciv_derived.ventilation`
  GROUP BY stay_id
),
cohort AS (
  SELECT
    b.*,
    s.sepsis_time,
    DATETIME_DIFF(vp.first_vp_time, s.sepsis_time, HOUR) AS first_vp_hours,
    CASE WHEN vp.first_vp_time >= s.sepsis_time
          AND vp.first_vp_time < DATETIME_ADD(s.sepsis_time, INTERVAL treatment_window_hours HOUR)
         THEN 1 ELSE 0 END AS early_vasopressor,
    sofa.sofa_24h,
    charlson.charlson,
    vitals.hr_mean, vitals.mbp_mean, vitals.rr_mean, vitals.temp_mean,
    chem.creatinine_mean, chem.lactate_mean, chem.wbc_mean, chem.platelet_mean,
    COALESCE(vent.mechanical_vent, 0) AS mechanical_vent
  FROM base b
  JOIN sepsis s ON b.stay_id = s.stay_id
  LEFT JOIN vp ON b.stay_id = vp.stay_id
  LEFT JOIN sofa ON b.stay_id = sofa.stay_id
  LEFT JOIN charlson ON b.subject_id = charlson.subject_id AND b.hadm_id = charlson.hadm_id
  LEFT JOIN vitals ON b.stay_id = vitals.stay_id
  LEFT JOIN chem ON b.stay_id = chem.stay_id
  LEFT JOIN vent ON b.stay_id = vent.stay_id
  WHERE b.age >= 18
    AND b.icu_stay_num = 1
    -- landmark: remain under observation at least treatment_window_hours after sepsis time
    AND b.outtime > DATETIME_ADD(s.sepsis_time, INTERVAL treatment_window_hours HOUR)
)
SELECT *
FROM cohort;
