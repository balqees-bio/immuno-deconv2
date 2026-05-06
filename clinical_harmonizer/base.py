"""
base.py — Abstract base class for all clinical harmonizers.

Subclasses implement only two methods:
  load(**kwargs)      → dict[str, pd.DataFrame]  (raw DataFrames by name)
  harmonize(raw)      → pd.DataFrame             (mapped to SCHEMA_COLUMNS)

Everything else — column ordering, validation, logging, CSV writing — is
inherited from ClinicalHarmonizer.
"""
import logging
from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# ── canonical output schema ───────────────────────────────────────────────────

SCHEMA_COLUMNS: list[str] = [
    # Dataset provenance
    "dataset_id",
    "dataset_name",
    "source_database",
    "accession_id",
    "publication_doi",
    "cohort_name",
    "disease_group",
    "cancer_type",
    "human_or_preclinical",
    # Patient identity & demographics
    "patient_id",
    "patient_id_original",
    "sex",
    "age_at_diagnosis",
    "vital_status",
    # Clinical features
    "primary_diagnosis",
    "cancer_subtype",
    "primary_site",
    "histology",
    "stage_overall",
    "grade",
    # Sample identity
    "sample_id",
    "sample_id_original",
    "sample_type",
    "specimen_type",
    "tumor_or_normal",
    "primary_or_metastatic",
    "collection_timepoint",
    # Treatment
    "treatment_received",
    "surgery_status",
    "chemotherapy_status",
    "radiotherapy_status",
    "immunotherapy_status",
    "targeted_therapy_status",
    # Overall survival
    "overall_survival_time",
    "overall_survival_unit",
    "overall_survival_event",
    # Progression-free / disease-free survival
    "progression_free_survival_time",
    "progression_free_survival_unit",
    "progression_free_survival_event",
    # Outcome flags
    "recurrence_status",
    "metastasis_status",
]


class ClinicalHarmonizer(ABC):
    """
    Base class for clinical metadata harmonizers.

    Usage
    -----
    harmonizer = SomeHarmonizer()
    df = harmonizer.run(
        output_path="results/clinical_harmonized.csv",
        patient_file="...",
        sample_file="...",
    )
    """

    # Subclasses may override to tag their source
    SOURCE_DATABASE: str = "unknown"

    # ── abstract interface ────────────────────────────────────────────────────

    @abstractmethod
    def load(self, **kwargs) -> dict[str, pd.DataFrame]:
        """
        Load raw source files and return them as named DataFrames.

        Returns
        -------
        dict with at minimum a 'patient' and/or 'sample' key
        """

    @abstractmethod
    def harmonize(self, raw: dict[str, pd.DataFrame]) -> pd.DataFrame:
        """
        Map raw fields → SCHEMA_COLUMNS.  All columns in SCHEMA_COLUMNS
        must be present in the returned DataFrame (fill with pd.NA when
        the source does not supply a field).
        """

    # ── public entry point ────────────────────────────────────────────────────

    def run(self, output_path: str | None = None, **load_kwargs) -> pd.DataFrame:
        """
        Load → harmonize → validate → (optionally) save.

        Parameters
        ----------
        output_path : destination CSV; skipped when None
        **load_kwargs: forwarded verbatim to self.load()

        Returns
        -------
        Harmonised DataFrame (samples × SCHEMA_COLUMNS)
        """
        logger.info("[Clinical] Loading raw files (%s)…", self.__class__.__name__)
        raw = self.load(**load_kwargs)

        logger.info("[Clinical] Harmonising to canonical schema…")
        df = self.harmonize(raw)

        # Ensure all schema columns exist (add missing ones as NA)
        for col in SCHEMA_COLUMNS:
            if col not in df.columns:
                df[col] = pd.NA

        df = df[SCHEMA_COLUMNS]         # enforce column order
        df = self._validate(df)

        if output_path is not None:
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, index=False)
            logger.info(
                "[Clinical] Saved %d samples × %d columns → %s",
                len(df), len(df.columns), output_path,
            )

        return df

    # ── validation ────────────────────────────────────────────────────────────

    @staticmethod
    def _validate(df: pd.DataFrame) -> pd.DataFrame:
        errors: list[str] = []

        # Required non-null columns
        for col in ("patient_id", "sample_id"):
            if col in df.columns:
                n_null = df[col].isna().sum()
                if n_null:
                    errors.append(f"'{col}' has {n_null} null values")

        # Survival events must be 0 or 1
        for col in ("overall_survival_event", "progression_free_survival_event"):
            if col in df.columns:
                bad = df[col].dropna()
                bad = bad[~bad.isin([0, 1])]
                if len(bad):
                    errors.append(
                        f"'{col}' contains non-0/1 values: {bad.unique().tolist()}"
                    )

        # Age sanity
        if "age_at_diagnosis" in df.columns:
            ages = pd.to_numeric(df["age_at_diagnosis"], errors="coerce").dropna()
            if len(ages) and not ages.between(0, 120).all():
                errors.append(
                    f"'age_at_diagnosis' has {(~ages.between(0, 120)).sum()} "
                    "values outside 0–120"
                )

        for e in errors:
            logger.error("[Clinical] Validation FAILED: %s", e)
        if not errors:
            logger.info("[Clinical] All validation checks passed.")

        # Log missing-rate per column
        missing = df.isna().mean().sort_values(ascending=False)
        high_missing = missing[missing > 0.5]
        if not high_missing.empty:
            logger.warning(
                "[Clinical] Columns with >50%% missing values:\n%s",
                high_missing.to_string(),
            )
        else:
            logger.info(
                "[Clinical] Max missing rate: %.1f%% (%s)",
                missing.iloc[0] * 100, missing.index[0],
            )

        return df
