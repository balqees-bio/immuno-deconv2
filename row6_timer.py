"""
row6_timer.py — Row 6: TIMER / TIMER2.0
Reference: Li et al. 2016 — 6 immune cell types calibrated on TCGA
Package: immunedeconv (timer method) via rpy2
Input: linear TPM matrix + cancer type label per sample

IMPORTANT: Cancer type strings must match TIMER's supported labels.
See immunedeconv::timer_available_cancers in R for the full list.
Using the wrong cancer type produces invalid results — scores are
calibrated against TCGA tumor-type-specific immune profiles.

Do NOT mix scores from the TIMER2.0 web portal with results from
this implementation — they use different background correction approaches.
"""
import logging
from pathlib import Path

import pandas as pd
import rpy2.robjects as ro

from rpy2_bridge import immunedeconv, df_to_r_matrix, r_to_df

logger = logging.getLogger(__name__)


class TIMERScorer:
    """Run TIMER deconvolution for one or multiple cancer types."""

    def run(
        self,
        expr: pd.DataFrame,
        output_dir: str,
        cancer_type: str | None = None,
        metadata_path: str | None = None,
    ) -> pd.DataFrame:
        """
        Parameters
        ----------
        expr          : linear TPM DataFrame, genes × samples, HGNC row index
        output_dir    : directory where CSV files are saved
        cancer_type   : TCGA cancer-type label (e.g. 'SKCM', 'BRCA', 'LUAD')
                        applied to ALL samples — use for single-cancer cohorts
        metadata_path : path to CSV with columns [sample_id, cancer_type]
                        — use for multi-cancer cohorts (overrides cancer_type)

        Returns
        -------
        DataFrame (6 cell types × samples, index = cell_type)
        """
        if cancer_type is None and metadata_path is None:
            raise ValueError("Provide either cancer_type (single) or metadata_path (multi-cancer).")

        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if metadata_path:
            return self._run_multi_cancer(expr, out, metadata_path)
        return self._run_single_cancer(expr, out, cancer_type)

    # ------------------------------------------------------------------
    def _run_single_cancer(
        self, expr: pd.DataFrame, out: Path, cancer_type: str
    ) -> pd.DataFrame:
        logger.info("Running TIMER for cancer type: %s (%d samples)...", cancer_type, expr.shape[1])

        r_mat = df_to_r_matrix(expr)
        r_result = immunedeconv.deconvolute(
            r_mat,
            "timer",
            indications=ro.StrVector([cancer_type] * expr.shape[1]),
        )
        timer = r_to_df(r_result).set_index("cell_type")
        fname = f"timer_{cancer_type}.csv"
        timer.to_csv(out / fname)

        logger.info(
            "  -> %d cell types × %d samples  [saved: %s]",
            timer.shape[0], timer.shape[1], fname,
        )
        self._validate(timer, cancer_type)
        return timer

    def _run_multi_cancer(
        self, expr: pd.DataFrame, out: Path, metadata_path: str
    ) -> pd.DataFrame:
        metadata = pd.read_csv(metadata_path, index_col=0)
        if "cancer_type" not in metadata.columns:
            raise ValueError("metadata CSV must contain a 'cancer_type' column.")

        results = []
        for ctype, samples in metadata.groupby("cancer_type"):
            shared = expr.columns.intersection(samples.index)
            if len(shared) == 0:
                logger.warning("No matching samples found for cancer type '%s' — skipping.", ctype)
                continue

            subset = expr[shared]
            logger.info(
                "Running TIMER for cancer type: %s (%d samples)...", ctype, len(shared)
            )
            r_mat = df_to_r_matrix(subset)
            r_result = immunedeconv.deconvolute(
                r_mat,
                "timer",
                indications=ro.StrVector([ctype] * len(shared)),
            )
            df = r_to_df(r_result).set_index("cell_type")
            results.append(df)
            logger.info("  -> %d cell types × %d samples", df.shape[0], df.shape[1])

        if not results:
            raise RuntimeError("TIMER produced no results — check metadata and sample IDs.")

        # Concatenate and restore original column order
        timer_all = pd.concat(results, axis=1)[expr.columns]
        timer_all.to_csv(out / "timer_all_cancertypes.csv")
        logger.info("Multi-cancer TIMER saved: timer_all_cancertypes.csv")
        self._validate(timer_all, "multi-cancer")
        return timer_all

    # ------------------------------------------------------------------
    @staticmethod
    def _validate(timer: pd.DataFrame, label: str) -> None:
        """Warn if any TIMER fractions are negative."""
        neg_count = int((timer < 0).sum().sum())
        if neg_count > 0:
            logger.warning(
                "TIMER (%s): %d negative fraction values detected. "
                "Verify cancer type label and input normalization.",
                label,
                neg_count,
            )
