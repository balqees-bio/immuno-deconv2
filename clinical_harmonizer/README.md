# Clinical Harmonizer

A source-agnostic framework for normalising clinical metadata from multiple databases into a single canonical schema, designed to accompany immune deconvolution and other multi-omics analyses.

---

## Contents

- [Overview](#overview)
- [Module Structure](#module-structure)
- [Canonical Schema](#canonical-schema)
- [Quick Start](#quick-start)
  - [Python API](#python-api)
  - [Pipeline CLI](#pipeline-cli)
- [TCGAHarmonizer](#tcgaharmonizer)
  - [Input File Format](#input-file-format)
  - [Raw → Schema Field Mapping](#raw--schema-field-mapping)
  - [Value Transformations](#value-transformations)
  - [Known Limitations (GDC Export)](#known-limitations-gdc-export)
- [Data Quality Report — TCGA BRCA](#data-quality-report--tcga-brca)
- [Pipeline Integration](#pipeline-integration)
  - [Row 0 — Clinical Harmonization](#row-0--clinical-harmonization)
  - [Row 7 — Automatic Clinical Join](#row-7--automatic-clinical-join)
  - [Output Files](#output-files)
- [Adding a New Source](#adding-a-new-source)
  - [Step-by-step guide](#step-by-step-guide)
  - [Example skeleton](#example-skeleton)
- [Validation Rules](#validation-rules)
- [Normalizer Reference](#normalizer-reference)

---

## Overview

Clinical data arrives from different databases (TCGA/GDC, cBioPortal, GEO, in-house LIMS) using different column names, value encodings, and file formats. This package provides:

- A **stable 41-column canonical schema** that all harmonized outputs conform to, regardless of source.
- An **abstract base class** (`ClinicalHarmonizer`) that enforces the schema contract, handles validation, missing-rate reporting, and CSV writing. Subclasses only implement two methods.
- A **`TCGAHarmonizer`** for GDC clinical exports, covering patient demographics, staging, survival endpoints, and sample metadata.
- **Automatic integration** with the Domain 2 deconvolution pipeline: once `clinical_harmonized.csv` exists in the results directory, `row7_merge.py` joins it onto the master deconvolution table automatically.

---

## Module Structure

```
clinical_harmonizer/
├── __init__.py   — public exports: ClinicalHarmonizer, TCGAHarmonizer, HARMONIZERS
├── base.py       — abstract base class + SCHEMA_COLUMNS definition + validation
├── tcga.py       — TCGAHarmonizer (GDC patient + sample files)
└── utils.py      — stateless Series-level normalizer functions
```

| File | Responsibility |
|------|---------------|
| `base.py` | Defines the schema, the abstract interface, and the shared `run()` / `_validate()` logic. Nothing source-specific lives here. |
| `tcga.py` | Reads GDC TSVs, handles the 4-row header format, maps all available raw fields, joins patient and sample tables. |
| `utils.py` | Pure functions: `normalize_sex`, `normalize_stage`, `survival_event_from_status`, etc. Reusable in any subclass. |
| `__init__.py` | Registers available harmonizers in `HARMONIZERS` dict so `main.py` can look them up by name. |

---

## Canonical Schema

Every harmonized output has exactly these 41 columns, in this order. Columns that a given source cannot populate are filled with `NaN`/`pd.NA`.

### Dataset Provenance

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `dataset_id` | str | Unique identifier for the dataset, prefixed by source | `GDC_TCGA-BRCA` |
| `dataset_name` | str | Human-readable dataset name | `Invasive Breast Carcinoma` |
| `source_database` | str | Database the data was downloaded from | `GDC` |
| `accession_id` | str | Database accession / project identifier | `TCGA-BRCA` |
| `publication_doi` | str | DOI of the primary publication; `NaN` if not set | `10.1038/nature11412` |
| `cohort_name` | str | Short cohort label, typically TCGA cancer type code | `BRCA` |
| `disease_group` | str | Broad disease category | `Cancer` |
| `cancer_type` | str | Cancer type code | `BRCA` |
| `human_or_preclinical` | str | `human` or `preclinical` | `human` |

### Patient Identity & Demographics

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `patient_id` | str | **Required.** Standardised patient identifier | `TCGA-3C-AAAU` |
| `patient_id_original` | str | Original/source patient identifier (e.g. GDC UUID) | `6e7d5ec6-…` |
| `sex` | str | `F` or `M` | `F` |
| `age_at_diagnosis` | float | Age in years at diagnosis | `55.0` |
| `vital_status` | str | `alive` or `dead` | `alive` |

### Clinical Features

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `primary_diagnosis` | str | WHO/ICD-O diagnosis text | `Infiltrating Lobular Carcinoma` |
| `cancer_subtype` | str | Oncotree or source-specific subtype | `Breast Invasive Lobular Carcinoma` |
| `primary_site` | str | Organ of origin | `Breast` |
| `histology` | str | ICD-O morphology code | `8520/3` |
| `stage_overall` | str | AJCC stage (prefix "Stage " stripped); `NaN` if indeterminate | `IIA` |
| `grade` | str | Tumour grade; `NaN` — not available in GDC export | |

### Sample Identity

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `sample_id` | str | **Required.** Standardised sample identifier; join key to expression matrix | `TCGA-3C-AAAU-01A` |
| `sample_id_original` | str | Original sample identifier | |
| `sample_type` | str | Raw sample type label from source | `Primary Tumor` |
| `specimen_type` | str | Oncotree code for the specimen | `BRCA` |
| `tumor_or_normal` | str | `tumor` or `normal` | `tumor` |
| `primary_or_metastatic` | str | `primary`, `metastatic`, or `recurrent` | `primary` |
| `collection_timepoint` | str | Timepoint relative to treatment; `NaN` — not in GDC export | |

### Treatment

| Column | Type | Description | Notes |
|--------|------|-------------|-------|
| `treatment_received` | str | `yes` or `no` — any prior treatment | GDC only provides aggregate flag |
| `surgery_status` | str | Surgery received | `NaN` — not in GDC export |
| `chemotherapy_status` | str | Chemotherapy received | `NaN` — not in GDC export |
| `radiotherapy_status` | str | Radiotherapy received | `NaN` — not in GDC export |
| `immunotherapy_status` | str | Immunotherapy received | `NaN` — not in GDC export |
| `targeted_therapy_status` | str | Targeted therapy received | `NaN` — not in GDC export |

### Overall Survival

| Column | Type | Allowed values | Notes |
|--------|------|----------------|-------|
| `overall_survival_time` | float | ≥ 0 | Duration in `overall_survival_unit` |
| `overall_survival_unit` | str | `months` | Fixed for GDC exports |
| `overall_survival_event` | Int64 | `0` (censored) / `1` (event: death) | Nullable integer |

### Progression-Free / Disease-Free Survival

| Column | Type | Allowed values | Notes |
|--------|------|----------------|-------|
| `progression_free_survival_time` | float | ≥ 0 | Duration in `progression_free_survival_unit` |
| `progression_free_survival_unit` | str | `months` | Fixed for GDC exports |
| `progression_free_survival_event` | Int64 | `0` (disease-free) / `1` (event: recurrence/progression) | Nullable integer |

### Outcome Flags

| Column | Type | Allowed values | Description |
|--------|------|----------------|-------------|
| `recurrence_status` | str | `recurred`, `disease_free` | Derived from DFS/PFS status field |
| `metastasis_status` | str | `metastatic`, `no_metastasis` | Derived from AJCC M-stage |

---

## Quick Start

### Python API

**Harmonize a TCGA GDC cohort:**

```python
from clinical_harmonizer import TCGAHarmonizer

harmonizer = TCGAHarmonizer()
df = harmonizer.run(
    output_path="results/clinical_harmonized.csv",
    patient_file="input-data/tcga_brca_gdc_data_clinical_patient.txt",
    sample_file="input-data/tcga_brca_gdc_data_clinical_sample.txt",
)
# df: 1102 rows × 41 columns
```

**Use without writing a file (in-memory only):**

```python
df = harmonizer.run(
    patient_file="...",
    sample_file="...",
    # output_path omitted → no file written
)
```

**Override dataset-level constants:**

```python
harmonizer = TCGAHarmonizer(
    publication_doi="10.1038/nature11412",
    human_or_preclinical="human",
    disease_group="Cancer",
)
```

**Access the schema column list programmatically:**

```python
from clinical_harmonizer import SCHEMA_COLUMNS
print(SCHEMA_COLUMNS)   # list of 41 column names in canonical order
```

### Pipeline CLI

Run clinical harmonization as part of the full Domain 2 pipeline:

```bash
python main.py \
  --input  expression_matrix.tsv \
  --output ./results \
  --cancer-type BRCA \
  --clinical-patient input-data/tcga_brca_gdc_data_clinical_patient.txt \
  --clinical-sample  input-data/tcga_brca_gdc_data_clinical_sample.txt \
  --clinical-source  tcga
```

**Run clinical harmonization alone** (skip all deconvolution rows):

```bash
python main.py \
  --input  expression_matrix.tsv \
  --output ./results \
  --clinical-patient input-data/tcga_brca_gdc_data_clinical_patient.txt \
  --clinical-sample  input-data/tcga_brca_gdc_data_clinical_sample.txt \
  --clinical-source  tcga \
  --skip-row3 --skip-row4 --skip-row5 --skip-row6
```

**Skip clinical harmonization** (deconvolution only):

```bash
python main.py \
  --input  expression_matrix.tsv \
  --output ./results \
  --cancer-type BRCA \
  --skip-clinical
```

---

## TCGAHarmonizer

### Input File Format

GDC clinical exports are tab-separated with **4 metadata header rows** before the data:

```
row 0  — human-readable column labels
row 1  — column descriptions
row 2  — data types (STRING / NUMBER)
row 3  — priority flags (0 or 1)
row 4  — machine-readable column names  ← used as DataFrame headers
row 5+ — patient/sample data
```

The harmonizer calls `pd.read_csv(..., skiprows=4)` to skip rows 0–3 and treat row 4 as the header. The following GDC sentinel values are all treated as `NaN`:

```
[Not Available]  [Not Applicable]  [Discrepancy]
[Unknown]        not reported      unknown
```

**Patient file** (`tcga_*_clinical_patient.txt`) — one row per patient:

```
PATIENT_ID  AGE  SEX  VITAL_STATUS  OS_STATUS  OS_MONTHS
DFS_STATUS  DFS_MONTHS  PATH_STAGE  PATH_M_STAGE
PRIMARY_DIAGNOSIS  PRIMARY_SITE_PATIENT  MORPHOLOGY
PRIOR_TREATMENT  PROJECT_ID  PROJECT_NAME  …
```

**Sample file** (`tcga_*_clinical_sample.txt`) — one row per sample (a patient may have multiple samples):

```
PATIENT_ID  SAMPLE_ID  SAMPLE_TYPE  ONCOTREE_CODE
CANCER_TYPE_DETAILED  OTHER_SAMPLE_ID  …
```

The two files are joined on `PATIENT_ID`. Every sample row gets all patient fields broadcast to it. Samples with no matching patient record are retained with patient columns set to `NaN` and a warning is logged.

---

### Raw → Schema Field Mapping

#### Patient file → schema

| Raw column | Schema column | Transform applied |
|------------|---------------|-------------------|
| `PROJECT_ID` | `dataset_id` | Prefixed with `"GDC_"` |
| `PROJECT_ID` | `accession_id` | As-is |
| `PROJECT_ID` | `cohort_name` | Split on `-`, take index `[1]` → `"BRCA"` |
| `PROJECT_ID` | `cancer_type` | Same as `cohort_name` |
| `PROJECT_NAME` | `dataset_name` | As-is |
| `OTHER_PATIENT_ID` | `patient_id_original` | GDC UUID |
| `PATIENT_ID` | `patient_id` | Strip whitespace |
| `SEX` | `sex` | `normalize_sex()` → `F` / `M` |
| `AGE` | `age_at_diagnosis` | Cast to float (already integer years in GDC export) |
| `VITAL_STATUS` | `vital_status` | `normalize_vital_status()` → lowercase |
| `PRIMARY_DIAGNOSIS` | `primary_diagnosis` | As-is |
| `PRIMARY_SITE_PATIENT` | `primary_site` | As-is |
| `MORPHOLOGY` | `histology` | ICD-O code, as-is |
| `PATH_STAGE` | `stage_overall` | `normalize_stage()` strips `"Stage "` prefix |
| `PATH_M_STAGE` | `metastasis_status` | `normalize_metastasis()` |
| `PRIOR_TREATMENT` | `treatment_received` | `normalize_prior_treatment()` → `yes` / `no` |
| `OS_MONTHS` | `overall_survival_time` | Cast to float |
| `OS_STATUS` | `overall_survival_event` | `survival_event_from_status()` → `0` / `1` |
| `DFS_MONTHS` | `progression_free_survival_time` | Cast to float |
| `DFS_STATUS` | `progression_free_survival_event` | `survival_event_from_status()` → `0` / `1` |
| `DFS_STATUS` | `recurrence_status` | `recurrence_from_dfs_status()` → `recurred` / `disease_free` |
| *(static)* | `source_database` | `"GDC"` |
| *(static)* | `overall_survival_unit` | `"months"` |
| *(static)* | `progression_free_survival_unit` | `"months"` |
| *(static)* | `human_or_preclinical` | `"human"` (constructor default) |
| *(static)* | `disease_group` | `"Cancer"` (constructor default) |
| *(unavailable)* | `grade`, `surgery_status`, `chemotherapy_status`, `radiotherapy_status`, `immunotherapy_status`, `targeted_therapy_status`, `publication_doi` | `NaN` |

#### Sample file → schema

| Raw column | Schema column | Transform applied |
|------------|---------------|-------------------|
| `PATIENT_ID` | `patient_id` | Join key (strip whitespace) |
| `SAMPLE_ID` | `sample_id` | Strip whitespace |
| `OTHER_SAMPLE_ID` | `sample_id_original` | As-is |
| `SAMPLE_TYPE` | `sample_type` | As-is (`"Primary Tumor"`, `"Metastatic"`) |
| `SAMPLE_TYPE` | `tumor_or_normal` | `tumor_or_normal()` regex → `"tumor"` / `"normal"` |
| `SAMPLE_TYPE` | `primary_or_metastatic` | `primary_or_metastatic()` → `"primary"` / `"metastatic"` |
| `ONCOTREE_CODE` | `specimen_type` | As-is |
| `CANCER_TYPE_DETAILED` | `cancer_subtype` | Oncotree detailed label |
| *(unavailable)* | `collection_timepoint` | `NaN` |

---

### Value Transformations

#### Sex

| Raw value | Output |
|-----------|--------|
| `Female` / `female` / `f` / `woman` | `F` |
| `Male` / `male` / `m` / `man` | `M` |
| anything else / NaN | `NaN` |

#### Stage

Strips the `"Stage "` prefix (case-insensitive). Values that indicate an indeterminate stage are mapped to `NaN`:

| Raw value | Output |
|-----------|--------|
| `Stage IIA` | `IIA` |
| `Stage IA` | `IA` |
| `Stage IV` | `IV` |
| `Stage X` | `NaN` |
| `[Not Available]` / `Unknown` | `NaN` |

#### Survival Events

GDC and cBioPortal use `"0:LIVING"` / `"1:DECEASED"` and `"0:DiseaseFree"` / `"1:Recurred/Progressed"` encodings. The parser handles:
- Strings matching `^\d+:` — extracts the leading digit as the event flag
- Plain `0` / `1` integers — used directly
- Anything else — `NaN` (nullable integer `pd.NA`)

| Raw OS_STATUS | `overall_survival_event` |
|---------------|--------------------------|
| `0:LIVING` | `0` |
| `1:DECEASED` | `1` |
| `NaN` | `<NA>` |

| Raw DFS_STATUS | `progression_free_survival_event` | `recurrence_status` |
|----------------|-----------------------------------|---------------------|
| `0:DiseaseFree` | `0` | `disease_free` |
| `1:Recurred/Progressed` | `1` | `recurred` |
| `NaN` | `<NA>` | `NaN` |

#### Metastasis (from M-stage)

| Raw PATH_M_STAGE | `metastasis_status` |
|------------------|---------------------|
| `M0` | `no_metastasis` |
| `M1` / `M1a` / `M1b` / `M1c` | `metastatic` |
| `MX` / `M0 (i+)` / `cM0 (I+)` | `NaN` |

#### Tumor / Normal (from sample type)

| `SAMPLE_TYPE` contains | `tumor_or_normal` | `primary_or_metastatic` |
|------------------------|-------------------|-------------------------|
| `tumor` / `metastat` / `cancer` | `tumor` | |
| `primary` / `initial` | | `primary` |
| `metastat` | | `metastatic` |
| `recur` / `relaps` | | `recurrent` |
| `normal` / `benign` / `blood` / `control` | `normal` | |

---

### Known Limitations (GDC Export)

| Field | Status | Reason |
|-------|--------|--------|
| `grade` | Not available | Not included in GDC standard clinical export TSVs |
| `surgery_status` | Not available | GDC only provides aggregate `PRIOR_TREATMENT` flag |
| `chemotherapy_status` | Not available | Same as above |
| `radiotherapy_status` | Not available | Same as above |
| `immunotherapy_status` | Not available | Same as above |
| `targeted_therapy_status` | Not available | Same as above |
| `collection_timepoint` | Not available | Not in GDC clinical TSVs |
| `publication_doi` | Optional | Pass via constructor: `TCGAHarmonizer(publication_doi="10.1038/...")` |
| DFS data | ~14% missing | Not all patients have disease-free survival follow-up in TCGA |

---

## Data Quality Report — TCGA BRCA

Produced from `tcga_brca_gdc_data_clinical_patient.txt` (1098 patients) and `tcga_brca_gdc_data_clinical_sample.txt` (1102 samples).

**Output:** 1102 rows × 41 columns

### Completeness by column

| Column | Completeness | Notes |
|--------|-------------|-------|
| `patient_id`, `sample_id`, `source_database`, `cancer_type`, `dataset_id` | **100%** | Always populated |
| `sex` | **99.9%** | 1 patient missing |
| `vital_status` | **99.9%** | 1 patient missing |
| `overall_survival_time` / `_event` | **99.9%** | 1 patient missing |
| `treatment_received` | **99.7%** | 3 patients missing |
| `primary_diagnosis`, `primary_site`, `histology` | **99.8–99.9%** | |
| `stage_overall` | **97.7%** | 25 samples — Stage X mapped to NaN |
| `age_at_diagnosis` | **98.5%** | 16 patients missing |
| `metastasis_status` | **84.6%** | MX coded as indeterminate → NaN |
| `progression_free_survival_time` / `_event` | **85.5–85.6%** | ~14% no DFS follow-up |
| `recurrence_status` | **85.6%** | Mirrors DFS availability |
| `grade`, `surgery_status`, `chemotherapy_status`, `radiotherapy_status`, `immunotherapy_status`, `targeted_therapy_status`, `collection_timepoint`, `publication_doi` | **0%** | Not available in GDC export — `NaN` by design |

### Value distributions (BRCA)

| Column | Values |
|--------|--------|
| `sex` | F: 1089 · M: 12 |
| `vital_status` | alive: 947 · dead: 154 |
| `stage_overall` | IIA: 358 · IIB: 260 · IIIA: 154 · I: 90 · IA: 86 · IIIC: 67 · IIIB: 27 · IV: 20 · IB: 7 |
| `tumor_or_normal` | tumor: 1102 (all BRCA samples are tumour) |
| `primary_or_metastatic` | primary: 1095 · metastatic: 7 |
| `recurrence_status` | disease_free: 861 · recurred: 82 |
| `metastasis_status` | no_metastasis: 910 · metastatic: 22 |
| `treatment_received` | no: 1086 · yes: 13 |

---

## Pipeline Integration

### Row 0 — Clinical Harmonization

Clinical harmonization runs as **Row 0** in `main.py`, before any expression analysis. It is triggered only when both `--clinical-patient` and `--clinical-sample` are supplied.

```
Row 0  → clinical_harmonized.csv
Row 3  → deconv_epic.csv, deconv_quantiseq.csv, deconv_mcp_counter.csv
Row 4  → xcell_cell_scores.csv, xcell_composite_scores.csv
Row 5  → ips_scores.csv
Row 6  → timer_{CANCER_TYPE}.csv
Row 7  → Domain2_master_scores.csv
         Domain2_master_with_clinical.csv  (auto, if clinical present)
```

### Row 7 — Automatic Clinical Join

`MasterMerger.run()` in `row7_merge.py` checks for `clinical_harmonized.csv` in the output directory automatically. If found, it left-joins on `sample_id` (clinical) = row index (deconvolution master) and writes a second output:

```
Domain2_master_with_clinical.csv   (samples × deconv features + clinical features)
```

The join is a **left join** keyed on `sample_id`, meaning:
- Deconvolution samples with no clinical match keep their deconv scores (clinical columns = NaN).
- Clinical samples with no matching expression data are excluded.
- Any clinical columns already present in the deconvolution master are dropped before joining to avoid collisions (logged at DEBUG level).

### Output Files

| File | Rows | Columns | Description |
|------|------|---------|-------------|
| `clinical_harmonized.csv` | 1 per sample | 41 | Canonical schema, all sources normalised |
| `Domain2_master_scores.csv` | 1 per sample | ~41 | Immune deconvolution features only |
| `Domain2_master_with_clinical.csv` | 1 per sample | ~81 | Deconvolution + clinical, joined on `sample_id` |

---

## Adding a New Source

To support a new clinical data source (e.g. cBioPortal, GEO, in-house LIMS), implement a subclass of `ClinicalHarmonizer`.

### Step-by-step guide

1. **Create** `clinical_harmonizer/my_source.py`
2. **Subclass** `ClinicalHarmonizer` and implement `load()` and `harmonize()`
3. **Use utils** from `utils.py` for all field normalizations — do not reimplement them
4. **Fill unavailable fields with `pd.NA`** — the base class inserts any missing schema columns automatically, but it is clearer to be explicit in `harmonize()`
5. **Register** the new class in `__init__.py` inside `HARMONIZERS`
6. **Test** with `harmonizer.run(output_path=None, ...)` and inspect the returned DataFrame

### Example skeleton

```python
# clinical_harmonizer/cbio.py
import pandas as pd
from .base import ClinicalHarmonizer
from .utils import normalize_sex, normalize_stage, survival_event_from_status

class CBioPortalHarmonizer(ClinicalHarmonizer):
    """Harmonizer for cBioPortal clinical data exports."""

    SOURCE_DATABASE = "cBioPortal"

    def load(self, clinical_file: str) -> dict[str, pd.DataFrame]:
        # cBioPortal uses a single combined file with a 4-row comment header
        df = pd.read_csv(clinical_file, sep="\t", comment="#", dtype=str)
        return {"combined": df}

    def harmonize(self, raw: dict) -> pd.DataFrame:
        df = raw["combined"].copy()
        out = pd.DataFrame()

        # --- dataset provenance ---
        out["source_database"]  = "cBioPortal"
        out["dataset_id"]       = "CBIO_" + df.get("STUDY_ID", pd.NA)
        out["human_or_preclinical"] = "human"
        # ... map remaining fields

        # --- patient ---
        out["patient_id"]       = df["PATIENT_ID"]
        out["sex"]              = normalize_sex(df["SEX"])
        out["age_at_diagnosis"] = pd.to_numeric(df.get("AGE"), errors="coerce")
        out["stage_overall"]    = normalize_stage(df.get("AJCC_PATHOLOGIC_TUMOR_STAGE", pd.NA))

        # --- survival ---
        out["overall_survival_time"]  = pd.to_numeric(df.get("OS_MONTHS"), errors="coerce")
        out["overall_survival_unit"]  = "months"
        out["overall_survival_event"] = survival_event_from_status(df.get("OS_STATUS", pd.NA), "1:")

        # --- fields not available ---
        out["grade"] = pd.NA
        # ... other unavailable fields

        return out
```

**Register in `__init__.py`:**

```python
from .cbio import CBioPortalHarmonizer

HARMONIZERS: dict[str, type[ClinicalHarmonizer]] = {
    "tcga":  TCGAHarmonizer,
    "cbio":  CBioPortalHarmonizer,   # ← add here
}
```

**Use from CLI:**

```bash
python main.py --clinical-source cbio --clinical-patient clinical_data.txt ...
```

---

## Validation Rules

`_validate()` runs automatically inside `run()` after harmonization. All checks log errors via the standard `logging` module (not raised as exceptions, so the pipeline continues).

| Rule | Severity | Check |
|------|----------|-------|
| `patient_id` never null | ERROR | `df["patient_id"].isna().sum() == 0` |
| `sample_id` never null | ERROR | `df["sample_id"].isna().sum() == 0` |
| Survival events are 0 or 1 | ERROR | `overall_survival_event` and `progression_free_survival_event` contain only `{0, 1, NaN}` |
| Age in physiological range | ERROR | `age_at_diagnosis` ∈ [0, 120] where not null |
| Missing rate per column | WARNING | Any column with > 50% missing values is listed |
| Max missing rate reported | INFO | Always printed so analysts can assess completeness |

---

## Normalizer Reference

All functions in `utils.py` accept a `pd.Series` and return a `pd.Series`. They handle `NaN` / `None` / `"nan"` input gracefully.

| Function | Input | Output |
|----------|-------|--------|
| `normalize_sex(s)` | Free-text sex labels | `"F"` / `"M"` / `NaN` |
| `normalize_stage(s)` | Stage strings with or without `"Stage "` prefix | Stage code (`"IIA"`, `"IV"`) / `NaN` for indeterminate |
| `normalize_vital_status(s)` | Free-text vital status | Lowercase string / `NaN` |
| `survival_event_from_status(s, dead_prefix)` | `"0:..."` / `"1:..."` or plain `0`/`1` | Nullable `Int64`: `0` / `1` / `<NA>` |
| `recurrence_from_dfs_status(s)` | DFS status string | `"recurred"` / `"disease_free"` / `NaN` |
| `normalize_metastasis(s)` | AJCC M-stage string | `"metastatic"` / `"no_metastasis"` / `NaN` |
| `tumor_or_normal(sample_type)` | Sample type label | `"tumor"` / `"normal"` / `NaN` |
| `primary_or_metastatic(sample_type)` | Sample type label | `"primary"` / `"metastatic"` / `"recurrent"` / `NaN` |
| `normalize_prior_treatment(s)` | `True`/`False`/`yes`/`no` | `"yes"` / `"no"` / `NaN` |
