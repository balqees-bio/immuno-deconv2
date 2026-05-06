"""
row7_merge.py — Row 7: Combine all Domain 2 outputs into a master table.

Loads the individual CSVs produced by Rows 3–6, transposes each so that
samples are rows, adds method-specific prefixes, and saves the merged
table as Domain2_master_scores.csv.

Column groups in the master table:
    EPIC_*      : 8 columns   — EPIC cell-type fractions (0–1)
    qT_*        : 11 columns  — quanTIseq fractions (0–1, sum ≈ 1)
    MCP_*       : 10 columns  — MCP-counter arbitrary units
    xCell_*     : 3 columns   — ImmuneScore / StromaScore / MicroenvironmentScore
    IPS         : 1 column    — total IPS score (0–10)
    IPS_class   : 1 column    — immunophenotype class A–D
    TIMER_*     : 6 columns   — TCGA-calibrated immune fractions
"""
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class MasterMerger:
    """Merge all Domain 2 CSVs into Domain2_master_scores.csv."""

    def run(
        self,
        output_dir: str,
        timer_cancer_type: str | None = None,
    ) -> pd.DataFrame:
        """
        Parameters
        ----------
        output_dir        : directory that contains all individual CSV outputs
        timer_cancer_type : cancer type string used for the single-cancer TIMER run
                            (e.g. 'SKCM'). Leave None if multi-cancer was used.

        Returns
        -------
        Master DataFrame (samples × ~40 features)
        """
        out = Path(output_dir)
        parts: list[pd.DataFrame] = []

        # ── Row 3 ─────────────────────────────────────────────────────────
        epic = self._load(out / "deconv_epic.csv", "EPIC")
        qt = self._load(out / "deconv_quantiseq.csv", "quanTIseq")
        mcp = self._load(out / "deconv_mcp_counter.csv", "MCP-counter")

        if epic is not None:
            parts.append(epic.T.add_prefix("EPIC_"))
        if qt is not None:
            parts.append(qt.T.add_prefix("qT_"))
        if mcp is not None:
            parts.append(mcp.T.add_prefix("MCP_"))

        # ── Row 4 — xCell composite scores only ───────────────────────────
        xcell_c = self._load(out / "xcell_composite_scores.csv", "xCell-composite")
        if xcell_c is not None:
            parts.append(xcell_c.T.add_prefix("xCell_"))

        # ── Row 5 — IPS ───────────────────────────────────────────────────
        ips = self._load_ips(out / "ips_scores.csv")
        if ips is not None:
            parts.append(ips)

        # ── Row 6 — TIMER ─────────────────────────────────────────────────
        timer_path = (
            out / f"timer_{timer_cancer_type}.csv"
            if timer_cancer_type
            else out / "timer_all_cancertypes.csv"
        )
        timer = self._load(timer_path, "TIMER")
        if timer is not None:
            parts.append(timer.T.add_prefix("TIMER_"))

        if not parts:
            raise RuntimeError("No output files found in %s — run the pipeline first." % out)

        master = pd.concat(parts, axis=1)
        master_path = out / "Domain2_master_scores.csv"
        master.to_csv(master_path)
        logger.info(
            "Master table saved: %s  (%d samples × %d features)",
            master_path,
            master.shape[0],
            master.shape[1],
        )

        self._validate(master)

        # ── Optional clinical join ─────────────────────────────────────────
        clinical_path = out / "clinical_harmonized.csv"
        if clinical_path.exists():
            master = self._join_clinical(master, clinical_path)

        return master

    def _join_clinical(self, master: pd.DataFrame, clinical_path) -> pd.DataFrame:
        """Left-join harmonised clinical metadata onto the master deconv table."""
        clinical = pd.read_csv(clinical_path, low_memory=False)

        # Align on sample_id (clinical) == master index (sample barcode)
        if "sample_id" not in clinical.columns:
            logger.warning("clinical_harmonized.csv has no 'sample_id' column — skipping join.")
            return master

        clinical_indexed = clinical.set_index("sample_id")

        # Drop columns already in master to avoid collisions
        overlap = master.columns.intersection(clinical_indexed.columns)
        if len(overlap):
            logger.debug("Dropping overlapping clinical columns before join: %s", overlap.tolist())
            clinical_indexed = clinical_indexed.drop(columns=overlap)

        combined = master.join(clinical_indexed, how="left")
        combined_path = clinical_path.parent / "Domain2_master_with_clinical.csv"
        combined.to_csv(combined_path)
        logger.info(
            "Master + clinical saved: %s  (%d samples × %d features)",
            combined_path, combined.shape[0], combined.shape[1],
        )
        return combined

    # ------------------------------------------------------------------
    def _load(self, path: Path, label: str) -> pd.DataFrame | None:
        """Load a cell_type × samples CSV; return None and warn if missing."""
        if not path.exists():
            logger.warning("Output file not found (skipping %s): %s", label, path)
            return None
        df = pd.read_csv(path, index_col=0)
        return df

    def _load_ips(self, path: Path) -> pd.DataFrame | None:
        """Load ips_scores.csv — already samples × columns format."""
        if not path.exists():
            logger.warning("Output file not found (skipping IPS): %s", path)
            return None
        ips = pd.read_csv(path, index_col=0)
        # Keep only the columns that go into the master table
        cols = [c for c in ["IPS_score", "Class"] if c in ips.columns]
        df = ips[cols].rename(columns={"IPS_score": "IPS", "Class": "IPS_class"})
        return df

    # ------------------------------------------------------------------
    @staticmethod
    def _validate(master: pd.DataFrame) -> None:
        """Run Section 7.3 validation checks and log results."""
        errors: list[str] = []

        if master.isnull().sum().sum() != 0:
            errors.append(
                f"NaN values in master table: {master.isnull().sum().sum()} total"
            )

        epic_cols = [c for c in master.columns if c.startswith("EPIC_")]
        if epic_cols:
            numeric_epic = master[epic_cols].select_dtypes(include="number")
            epic_sums = numeric_epic.sum(axis=1)
            over = (epic_sums > 1.05).sum()
            if over:
                errors.append(
                    f"EPIC fractions exceed 1.05 in {over} samples (max={epic_sums.max():.3f})"
                )

        if "IPS" in master.columns:
            ips_bad = (~master["IPS"].between(0, 10)).sum()
            if ips_bad:
                errors.append(f"IPS scores outside 0–10 range: {ips_bad} samples")

        timer_cols = [c for c in master.columns if c.startswith("TIMER_")]
        if timer_cols:
            numeric_timer = master[timer_cols].select_dtypes(include="number")
            neg = int((numeric_timer < 0).sum().sum())
            if neg:
                errors.append(f"TIMER has {neg} negative fraction values")

        xcell_me = "xCell_microenvironment score"
        if xcell_me in master.columns:
            neg_me = int((master[xcell_me] < 0).sum())
            if neg_me:
                errors.append(f"xCell MicroenvironmentScore < 0 in {neg_me} samples")

        if errors:
            for e in errors:
                logger.error("Validation FAILED: %s", e)
        else:
            logger.info("All validation checks passed.")
