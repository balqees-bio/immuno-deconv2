"""
Stateless normalizer functions shared across all harmonizer subclasses.
Each function accepts a raw value (or pandas Series) and returns a
standardised value (or Series).  All accept None/NaN gracefully.
"""
import re
import numpy as np
import pandas as pd

# ── sex ──────────────────────────────────────────────────────────────────────

_SEX_MAP = {
    "female": "F", "f": "F", "woman": "F",
    "male":   "M", "m": "M", "man":   "M",
}

def normalize_sex(s) -> pd.Series:
    return (
        s.astype(str).str.strip().str.lower()
         .map(_SEX_MAP)
    )


# ── stage ────────────────────────────────────────────────────────────────────

_STAGE_STRIP = re.compile(r"^stage\s+", re.IGNORECASE)
_STAGE_UNKNOWN = {"x", "unknown", "not reported", "nan", "none", ""}

def normalize_stage(s) -> pd.Series:
    cleaned = s.astype(str).str.strip()
    cleaned = cleaned.str.replace(_STAGE_STRIP, "", regex=True)
    return cleaned.where(~cleaned.str.lower().isin(_STAGE_UNKNOWN), other=pd.NA)


# ── vital status ─────────────────────────────────────────────────────────────

def normalize_vital_status(s) -> pd.Series:
    return s.astype(str).str.strip().str.lower().replace({"nan": pd.NA})


# ── survival event (integer 0/1) ─────────────────────────────────────────────

def survival_event_from_status(s: pd.Series, dead_prefix: str) -> pd.Series:
    """
    Convert cBioPortal/GDC status strings → 0/1 integer event flag.

    dead_prefix : the leading token that means event occurred,
                  e.g. '1:DECEASED' or '1:Recurred/Progressed'
                  (matched by checking whether the raw value starts with '1:')
    """
    numeric = pd.to_numeric(s, errors="coerce")
    valid_numeric = numeric.isin([0.0, 1.0])

    # Try to parse "0:..." / "1:..." pattern
    str_parsed = s.astype(str).str.extract(r"^(\d+):", expand=False)
    str_numeric = pd.to_numeric(str_parsed, errors="coerce")

    result = numeric.where(valid_numeric, str_numeric)
    return result.astype("Int64")   # nullable integer


# ── recurrence label ─────────────────────────────────────────────────────────

def recurrence_from_dfs_status(s: pd.Series) -> pd.Series:
    event = survival_event_from_status(s, "1:")
    labels = event.map({1: "recurred", 0: "disease_free"})
    return labels.where(event.notna(), other=pd.NA)


# ── metastasis from M-stage ───────────────────────────────────────────────────

_M_MAP = {
    "m0": "no_metastasis",
    "m1": "metastatic",
    "m1a": "metastatic", "m1b": "metastatic", "m1c": "metastatic",
}
_M_UNKNOWN = {"mx", "m0 (i+)", "cm0 (i+)", "nan", "none", "not reported", ""}

def normalize_metastasis(s: pd.Series) -> pd.Series:
    cleaned = s.astype(str).str.strip().str.lower()
    result = cleaned.map(_M_MAP)
    unknown_mask = cleaned.isin(_M_UNKNOWN)
    return result.where(~unknown_mask, other=pd.NA)


# ── sample type helpers ───────────────────────────────────────────────────────

def tumor_or_normal(sample_type: pd.Series) -> pd.Series:
    lower = sample_type.astype(str).str.lower()
    cond_tumor  = lower.str.contains("tumor|metastat|recur|cancer", na=False)
    cond_normal = lower.str.contains("normal|benign|blood|control", na=False)
    result = pd.Series(pd.NA, index=sample_type.index, dtype=object)
    result[cond_tumor]  = "tumor"
    result[cond_normal] = "normal"
    return result

def primary_or_metastatic(sample_type: pd.Series) -> pd.Series:
    lower = sample_type.astype(str).str.lower()
    result = pd.Series(pd.NA, index=sample_type.index, dtype=object)
    result[lower.str.contains("primary|initial", na=False)] = "primary"
    result[lower.str.contains("metastat", na=False)] = "metastatic"
    result[lower.str.contains("recur|relaps", na=False)] = "recurrent"
    return result


# ── prior treatment ───────────────────────────────────────────────────────────

def normalize_prior_treatment(s: pd.Series) -> pd.Series:
    lower = s.astype(str).str.strip().str.lower()
    result = pd.Series(pd.NA, index=s.index, dtype=object)
    result[lower.isin(["true",  "yes", "1"])] = "yes"
    result[lower.isin(["false", "no",  "0"])] = "no"
    return result
