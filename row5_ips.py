"""
row5_ips.py — Row 5: IPS (Immunophenogram Score)
Reference: Charoentong et al. 2017 Cell Reports — 122 genes, 4 modules
Package: immunedeconv (ips_rand method) via rpy2

IMPORTANT: IPS requires log2(TPM+1) input. This module accepts LINEAR TPM
and applies the log2 transformation internally to prevent the common mistake
of passing already-log-transformed data (which produces all Class D results).
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from rpy2_bridge import immunedeconv, df_to_r_matrix, r_to_df

logger = logging.getLogger(__name__)

_MODULES = ["MHC", "EC", "CP", "SC"]
_IPS_ROW = "IPS"


def _classify_ips(score: float) -> str:
    """Assign immunophenotype class A–D based on IPS total score."""
    if score >= 7:
        return "A"
    elif score >= 5:
        return "B"
    elif score >= 3:
        return "C"
    return "D"


class IPSScorer:
    """Score samples with the Immunophenogram (IPS) and assign class A–D."""

    def run(self, expr: pd.DataFrame, output_dir: str) -> pd.DataFrame:
        """
        Parameters
        ----------
        expr       : LINEAR TPM DataFrame (genes × samples, HGNC row index).
                     The log2 transform is applied here — do not pre-transform.
        output_dir : directory where ips_scores.csv is saved

        Returns
        -------
        DataFrame with columns: IPS_score, MHC, EC, CP, SC, Class
        (index = sample IDs)
        """
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        # IPS requires log2(TPM+1) — transform internally
        expr_log2 = np.log2(expr + 1)
        logger.info(
            "Running IPS (log2-transformed input: max=%.2f)...",
            float(expr_log2.values.max()),
        )

        r_mat = df_to_r_matrix(expr_log2)
        r_result = immunedeconv.deconvolute(r_mat, "ips_rand")
        ips_full = r_to_df(r_result).set_index("cell_type")

        # Extract total IPS and 4 module scores
        if _IPS_ROW not in ips_full.index:
            raise ValueError(
                f"'IPS' row not found in immunedeconv output. "
                f"Available rows: {ips_full.index.tolist()}"
            )

        ips_score = ips_full.loc[_IPS_ROW]
        module_scores = {m: ips_full.loc[m] for m in _MODULES if m in ips_full.index}
        missing_modules = [m for m in _MODULES if m not in ips_full.index]
        if missing_modules:
            logger.warning("IPS modules not found in output: %s", missing_modules)

        ips_class = ips_score.apply(_classify_ips)

        output = pd.DataFrame({"IPS_score": ips_score, **module_scores, "Class": ips_class})
        output.to_csv(out / "ips_scores.csv")

        class_counts = ips_class.value_counts().to_dict()
        logger.info(
            "  -> IPS scored for %d samples. Classes: %s",
            len(output),
            class_counts,
        )

        # Warn if scores are out of expected range
        ips_min = float(ips_score.min())
        ips_max = float(ips_score.max())
        if ips_min < 0 or ips_max > 10:
            logger.warning(
                "IPS scores outside 0–10 range (min=%.2f, max=%.2f). "
                "Verify input normalization.",
                ips_min,
                ips_max,
            )

        return output
