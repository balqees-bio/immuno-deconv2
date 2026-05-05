"""
row4_xcell.py — Row 4: xCell (64 cell types + composite TME scores)
Package: immunedeconv (xCell method) via rpy2
Input: linear TPM matrix (same matrix as Row 3 — do NOT log-transform)

xCell enrichment scores are NOT fractions and do NOT sum to 1.
The three composite scores (ImmuneScore, StromaScore, MicroenvironmentScore)
are xCell-exclusive and are split out into a separate CSV.
"""
import logging
from pathlib import Path

import pandas as pd

from rpy2_bridge import immunedeconv, df_to_r_matrix, r_to_df

logger = logging.getLogger(__name__)

_COMPOSITE_NAMES = {"ImmuneScore", "StromaScore", "MicroenvironmentScore"}


class XCellScorer:
    """Run xCell deconvolution and split results into cell scores + composite scores."""

    def run(self, expr: pd.DataFrame, output_dir: str) -> dict[str, pd.DataFrame]:
        """
        Parameters
        ----------
        expr       : linear TPM DataFrame, genes × samples, HGNC row index
        output_dir : directory where CSV files are saved

        Returns
        -------
        {
            'cell_scores' : DataFrame  (up to 64 cell types × samples),
            'composite'   : DataFrame  (3 composite scores × samples),
        }
        """
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        logger.info("Running xCell (64 cell types + composite scores)...")
        r_mat = df_to_r_matrix(expr)
        r_result = immunedeconv.deconvolute(r_mat, "xcell")
        xcell = r_to_df(r_result).set_index("cell_type")

        composite = xcell.loc[xcell.index.isin(_COMPOSITE_NAMES)]
        cell_scores = xcell.loc[~xcell.index.isin(_COMPOSITE_NAMES)]

        cell_scores.to_csv(out / "xcell_cell_scores.csv")
        composite.to_csv(out / "xcell_composite_scores.csv")

        logger.info(
            "  -> %d cell types scored; composite scores: %s",
            len(cell_scores),
            composite.index.tolist(),
        )

        # Sanity check: MicroenvironmentScore should be non-negative
        if "MicroenvironmentScore" in composite.index:
            me_min = float(composite.loc["MicroenvironmentScore"].min())
            if me_min < 0:
                logger.warning(
                    "MicroenvironmentScore contains negative values (min=%.4f). "
                    "Verify input is linear TPM (not log-transformed).",
                    me_min,
                )

        return {"cell_scores": cell_scores, "composite": composite}
