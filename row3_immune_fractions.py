"""
row3_immune_fractions.py — Row 3: Immune Cell Fractions
Tools: EPIC · quanTIseq · MCP-counter · CIBERSORTx (optional, token required)
Package: immunedeconv (Bioconductor) via rpy2
Input: linear TPM matrix (genes × samples, HGNC symbols)
"""
import logging
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

from rpy2_bridge import immunedeconv, df_to_r_matrix, r_to_df

logger = logging.getLogger(__name__)

# Cell-type name used by each method for the CD8+ T cell (for concordance check)
_CD8_LABEL = {
    "epic": "CD8+ T cell",
    "quantiseq": "T cell CD8+",
    "mcp_counter": "T cells.CD8",
}


class ImmuneFrections:
    """Run EPIC, quanTIseq, MCP-counter (and optionally CIBERSORTx) via immunedeconv."""

    def run(
        self,
        expr: pd.DataFrame,
        output_dir: str,
        methods: list[str] | None = None,
        cibersort_binary: str | None = None,
        cibersort_mat: str | None = None,
    ) -> dict[str, pd.DataFrame]:
        """
        Parameters
        ----------
        expr             : linear TPM DataFrame, genes × samples, HGNC row index
        output_dir       : directory where CSV files are saved
        methods          : list of immunedeconv method names to run
                           (default: ['epic', 'quantiseq', 'mcp_counter'])
        cibersort_binary : path to CIBERSORT.R (required for 'cibersort' / 'cibersort_abs')
        cibersort_mat    : path to LM22.txt signature matrix

        Returns
        -------
        dict mapping method name → DataFrame (cell_type × samples, index = cell_type)
        """
        if methods is None:
            methods = ["epic", "quantiseq", "mcp_counter"]

        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        # Configure CIBERSORTx if credentials supplied
        if cibersort_binary and cibersort_mat:
            logger.info("Configuring CIBERSORTx binary and signature matrix.")
            immunedeconv.set_cibersort_binary(cibersort_binary)
            immunedeconv.set_cibersort_mat(cibersort_mat)

        r_mat = df_to_r_matrix(expr)
        results: dict[str, pd.DataFrame] = {}

        for method in methods:
            logger.info("Running %s ...", method)
            try:
                r_result = immunedeconv.deconvolute(r_mat, method)
                df = r_to_df(r_result).set_index("cell_type")
                results[method] = df
                fname = f"deconv_{method}.csv"
                df.to_csv(out / fname)
                logger.info("  -> %d cell types × %d samples  [saved: %s]",
                            df.shape[0], df.shape[1], fname)
            except Exception as exc:
                logger.error("  %s failed: %s", method, exc)

        self._concordance_check(results)
        return results

    # ------------------------------------------------------------------
    def _concordance_check(self, results: dict[str, pd.DataFrame]) -> None:
        """Log Spearman rho between EPIC and quanTIseq CD8+ T cell estimates."""
        if "epic" not in results or "quantiseq" not in results:
            return

        cd8_epic_label = _CD8_LABEL["epic"]
        cd8_qt_label = _CD8_LABEL["quantiseq"]

        epic_df = results["epic"]
        qt_df = results["quantiseq"]

        if cd8_epic_label not in epic_df.index:
            logger.warning("CD8 concordance check skipped: '%s' not in EPIC output.", cd8_epic_label)
            return
        if cd8_qt_label not in qt_df.index:
            logger.warning("CD8 concordance check skipped: '%s' not in quanTIseq output.", cd8_qt_label)
            return

        shared = epic_df.columns.intersection(qt_df.columns)
        if len(shared) < 3:
            logger.warning("CD8 concordance check skipped: fewer than 3 shared samples.")
            return

        rho, pval = spearmanr(
            epic_df.loc[cd8_epic_label, shared].values,
            qt_df.loc[cd8_qt_label, shared].values,
        )
        level = logging.INFO if rho > 0.6 else logging.WARNING
        logger.log(level, "CD8 T concordance — EPIC vs quanTIseq: rho=%.3f  p=%.3e", rho, pval)
        if rho <= 0.6:
            logger.warning(
                "rho=%.3f is below the expected threshold of 0.6. "
                "Check input normalization or sample composition.",
                rho,
            )
