"""
tcga.py — TCGAHarmonizer: maps GDC patient + sample clinical exports
to the canonical SCHEMA_COLUMNS defined in base.py.

GDC clinical TSVs have 4 header rows before the data:
  row 0 — human-readable column labels
  row 1 — descriptions
  row 2 — data types
  row 3 — priority flags
  row 4 — machine column names  ← used as column headers (skiprows=4)
  row 5+ — data
"""
import logging

import pandas as pd

from .base import ClinicalHarmonizer
from .utils import (
    normalize_sex,
    normalize_stage,
    normalize_vital_status,
    survival_event_from_status,
    recurrence_from_dfs_status,
    normalize_metastasis,
    tumor_or_normal,
    primary_or_metastatic,
    normalize_prior_treatment,
)

logger = logging.getLogger(__name__)

# GDC exports use "not reported" / "[Not Available]" / "[Not Applicable]" as NA
_GDC_NA = {
    "[not available]", "[not applicable]", "[discrepancy]",
    "[unknown]", "not reported", "unknown", "nan", "",
}


def _read_gdc_tsv(path: str) -> pd.DataFrame:
    """Read a GDC TSV, using row 5 (index 4) as the header."""
    df = pd.read_csv(path, sep="\t", skiprows=4, dtype=str, na_values=list(_GDC_NA))
    df.columns = df.columns.str.strip()
    return df


class TCGAHarmonizer(ClinicalHarmonizer):
    """
    Harmonizer for TCGA / GDC clinical exports.

    Inputs
    ------
    patient_file : path to  tcga_*_clinical_patient.txt
    sample_file  : path to  tcga_*_clinical_sample.txt

    Dataset-level constants (auto-derived from PROJECT_ID / PROJECT_NAME
    in the patient file, or overridable via constructor).
    """

    SOURCE_DATABASE = "GDC"

    def __init__(
        self,
        publication_doi: str | None = None,
        human_or_preclinical: str = "human",
        disease_group: str = "Cancer",
    ):
        self._pub_doi = publication_doi
        self._human_or_preclinical = human_or_preclinical
        self._disease_group = disease_group

    # ── load ─────────────────────────────────────────────────────────────────

    def load(self, patient_file: str, sample_file: str) -> dict[str, pd.DataFrame]:
        pt = _read_gdc_tsv(patient_file)
        sp = _read_gdc_tsv(sample_file)
        logger.info(
            "[TCGA] Loaded patient=%s (%d rows), sample=%s (%d rows)",
            patient_file, len(pt), sample_file, len(sp),
        )
        return {"patient": pt, "sample": sp}

    # ── harmonize ────────────────────────────────────────────────────────────

    def harmonize(self, raw: dict[str, pd.DataFrame]) -> pd.DataFrame:
        pt = raw["patient"].copy()
        sp = raw["sample"].copy()

        # ── Patient table ─────────────────────────────────────────────────
        pt_out = pd.DataFrame(index=pt.index)

        # Dataset provenance (derived from project fields)
        project_id   = pt["PROJECT_ID"].fillna("UNKNOWN")
        project_name = pt.get("PROJECT_NAME", pd.Series("", index=pt.index)).fillna("")

        pt_out["dataset_id"]           = "GDC_" + project_id
        pt_out["dataset_name"]         = project_name
        pt_out["source_database"]      = "GDC"
        pt_out["accession_id"]         = project_id
        pt_out["publication_doi"]      = self._pub_doi
        pt_out["cohort_name"]          = project_id.str.split("-").str[1]
        pt_out["disease_group"]        = self._disease_group
        pt_out["cancer_type"]          = project_id.str.split("-").str[1]
        pt_out["human_or_preclinical"] = self._human_or_preclinical

        # Patient identity
        pt_out["patient_id"]          = pt["PATIENT_ID"].str.strip()
        pt_out["patient_id_original"] = pt.get("OTHER_PATIENT_ID", pd.NA)

        # Demographics
        pt_out["sex"]               = normalize_sex(pt["SEX"])
        pt_out["age_at_diagnosis"]  = pd.to_numeric(pt["AGE"], errors="coerce")
        pt_out["vital_status"]      = normalize_vital_status(pt["VITAL_STATUS"])

        # Clinical features
        pt_out["primary_diagnosis"] = pt.get("PRIMARY_DIAGNOSIS", pd.NA)
        pt_out["primary_site"]      = pt.get("PRIMARY_SITE_PATIENT", pd.NA)
        pt_out["histology"]         = pt.get("MORPHOLOGY", pd.NA)
        pt_out["stage_overall"]     = normalize_stage(
            pt.get("PATH_STAGE", pd.Series(pd.NA, index=pt.index))
        )
        pt_out["grade"] = pd.NA   # not in GDC export

        # Treatment (only aggregate prior-treatment flag in GDC)
        pt_out["treatment_received"]      = normalize_prior_treatment(
            pt.get("PRIOR_TREATMENT", pd.Series(pd.NA, index=pt.index))
        )
        pt_out["surgery_status"]          = pd.NA
        pt_out["chemotherapy_status"]     = pd.NA
        pt_out["radiotherapy_status"]     = pd.NA
        pt_out["immunotherapy_status"]    = pd.NA
        pt_out["targeted_therapy_status"] = pd.NA

        # Overall survival
        pt_out["overall_survival_time"]  = pd.to_numeric(pt.get("OS_MONTHS",  pd.NA), errors="coerce")
        pt_out["overall_survival_unit"]  = "months"
        pt_out["overall_survival_event"] = survival_event_from_status(
            pt.get("OS_STATUS", pd.Series(pd.NA, index=pt.index)), "1:"
        )

        # Progression-free / disease-free survival
        pt_out["progression_free_survival_time"]  = pd.to_numeric(pt.get("DFS_MONTHS", pd.NA), errors="coerce")
        pt_out["progression_free_survival_unit"]  = "months"
        pt_out["progression_free_survival_event"] = survival_event_from_status(
            pt.get("DFS_STATUS", pd.Series(pd.NA, index=pt.index)), "1:"
        )

        # Outcome flags
        pt_out["recurrence_status"] = recurrence_from_dfs_status(
            pt.get("DFS_STATUS", pd.Series(pd.NA, index=pt.index))
        )
        pt_out["metastasis_status"] = normalize_metastasis(
            pt.get("PATH_M_STAGE", pd.Series(pd.NA, index=pt.index))
        )

        # ── Sample table ──────────────────────────────────────────────────
        sp_out = pd.DataFrame(index=sp.index)
        sp_out["patient_id"]         = sp["PATIENT_ID"].str.strip()
        sp_out["sample_id"]          = sp["SAMPLE_ID"].str.strip()
        sp_out["sample_id_original"] = sp.get("OTHER_SAMPLE_ID", pd.NA)
        sp_out["sample_type"]        = sp.get("SAMPLE_TYPE", pd.NA)
        sp_out["specimen_type"]      = sp.get("ONCOTREE_CODE", pd.NA)
        sp_out["cancer_subtype"]     = sp.get("CANCER_TYPE_DETAILED", pd.NA)
        sp_out["tumor_or_normal"]    = tumor_or_normal(
            sp.get("SAMPLE_TYPE", pd.Series("", index=sp.index))
        )
        sp_out["primary_or_metastatic"] = primary_or_metastatic(
            sp.get("SAMPLE_TYPE", pd.Series("", index=sp.index))
        )
        sp_out["collection_timepoint"] = pd.NA

        # ── Join: one row per sample, patient fields broadcast ────────────
        merged = sp_out.merge(pt_out, on="patient_id", how="left")

        n_unmatched = merged["vital_status"].isna().sum()
        if n_unmatched:
            logger.warning(
                "[TCGA] %d samples had no matching patient record after join.",
                n_unmatched,
            )

        logger.info(
            "[TCGA] Harmonised: %d samples from %d patients",
            len(merged),
            merged["patient_id"].nunique(),
        )
        return merged
