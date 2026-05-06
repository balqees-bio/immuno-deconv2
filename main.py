"""
main.py — Domain 2 Immune Deconvolution Pipeline orchestrator.

Usage
-----
python main.py \\
    --input  expression_tpm.csv \\
    --output ./results \\
    --cancer-type SKCM \\
    [--clinical-patient clinical_patient.txt] \\
    [--clinical-sample  clinical_sample.txt] \\
    [--clinical-source  tcga] \\
    [--metadata metadata.csv] \\
    [--methods epic quantiseq mcp_counter] \\
    [--cibersort-binary /path/CIBERSORT.R] \\
    [--cibersort-mat /path/LM22.txt] \\
    [--skip-clinical] \\
    [--skip-row3] [--skip-row4] [--skip-row5] [--skip-row6]
"""
import argparse
import logging
import sys
from pathlib import Path


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s  %(levelname)-8s  %(name)s — %(message)s",
        datefmt="%H:%M:%S",
        level=level,
        stream=sys.stdout,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Domain 2 — Immune Deconvolution & Microenvironment Profiling"
    )
    p.add_argument("--input", required=True, metavar="FILE",
                   help="Expression matrix (CSV or TSV). Auto-detects orientation, "
                        "Ensembl IDs, and log2 encoding from the filename.")
    p.add_argument("--log2-offset", type=float, default=None, metavar="FLOAT",
                   help="Pseudocount used before log2 transform (e.g. 0.001 for "
                        "log2(TPM+0.001)). Auto-parsed from filename when possible; "
                        "use this flag to override.")
    p.add_argument("--output", default="./results", metavar="DIR",
                   help="Output directory (default: ./results)")
    p.add_argument("--cancer-type", metavar="TYPE",
                   help="TCGA cancer type label for TIMER (e.g. SKCM, BRCA, LUAD)")
    p.add_argument("--metadata", metavar="CSV",
                   help="Sample metadata CSV with 'cancer_type' column (for multi-cancer TIMER)")
    p.add_argument("--methods", nargs="+",
                   default=["epic", "quantiseq", "mcp_counter"],
                   metavar="METHOD",
                   help="Row 3 immunedeconv methods (default: epic quantiseq mcp_counter)")
    p.add_argument("--cibersort-binary", metavar="PATH",
                   help="Path to CIBERSORT.R (enables CIBERSORTx; also set --cibersort-mat)")
    p.add_argument("--cibersort-mat", metavar="PATH",
                   help="Path to LM22.txt signature matrix (required with --cibersort-binary)")
    # ── clinical harmonization ────────────────────────────────────────────
    p.add_argument("--clinical-patient", metavar="FILE",
                   help="Patient-level clinical file (e.g. GDC clinical_patient.txt)")
    p.add_argument("--clinical-sample", metavar="FILE",
                   help="Sample-level clinical file (e.g. GDC clinical_sample.txt)")
    p.add_argument("--clinical-source", default="tcga",
                   choices=["tcga"],
                   metavar="SOURCE",
                   help="Clinical data source / harmonizer to use (default: tcga)")
    p.add_argument("--skip-clinical", action="store_true",
                   help="Skip clinical harmonization (Row 0)")
    # ── visualization ─────────────────────────────────────────────────────
    p.add_argument("--skip-plots", action="store_true",
                   help="Skip Row 8 visualization (plots still runnable via plots.py)")
    p.add_argument("--dpi", type=int, default=150,
                   help="Plot resolution in DPI (default: 150)")
    # ── deconvolution rows ────────────────────────────────────────────────
    p.add_argument("--skip-row3", action="store_true",
                   help="Skip Row 3 (immune cell fractions)")
    p.add_argument("--skip-row4", action="store_true",
                   help="Skip Row 4 (xCell)")
    p.add_argument("--skip-row5", action="store_true",
                   help="Skip Row 5 (IPS)")
    p.add_argument("--skip-row6", action="store_true",
                   help="Skip Row 6 (TIMER)")
    p.add_argument("--verbose", action="store_true",
                   help="Enable DEBUG-level logging")
    return p.parse_args()


class Domain2Pipeline:
    """Orchestrates the four immune deconvolution steps for Domain 2."""

    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(self.__class__.__name__)

    def run(self) -> None:
        args = self.args

        # ── Row 0: Clinical harmonization ──────────────────────────────────
        if not args.skip_clinical:
            self._run_clinical()
        else:
            self.logger.info("Row 0 skipped (--skip-clinical).")

        # ── Load & validate expression matrix (only when deconv rows will run) ──
        _deconv_needed = not (args.skip_row3 and args.skip_row4
                              and args.skip_row5 and args.skip_row6)
        expr_tpm = None
        if _deconv_needed:
            self.logger.info("Loading expression matrix: %s", args.input)
            from input_validation import load_expression_matrix
            expr_tpm, _ = load_expression_matrix(args.input, log2_offset=args.log2_offset)

        # ── Row 3: Immune Cell Fractions ───────────────────────────────────
        if not args.skip_row3:
            self.logger.info("=" * 60)
            self.logger.info("Row 3 — Immune Cell Fractions")
            from row3_immune_fractions import ImmuneFrections
            ImmuneFrections().run(
                expr=expr_tpm,
                output_dir=str(self.output_dir),
                methods=args.methods,
                cibersort_binary=args.cibersort_binary,
                cibersort_mat=args.cibersort_mat,
            )
        else:
            self.logger.info("Row 3 skipped (--skip-row3).")

        # ── Row 4: xCell ───────────────────────────────────────────────────
        if not args.skip_row4:
            self.logger.info("=" * 60)
            self.logger.info("Row 4 — xCell")
            from row4_xcell import XCellScorer
            XCellScorer().run(expr=expr_tpm, output_dir=str(self.output_dir))
        else:
            self.logger.info("Row 4 skipped (--skip-row4).")

        # ── Row 5: IPS ─────────────────────────────────────────────────────
        if not args.skip_row5:
            self.logger.info("=" * 60)
            self.logger.info("Row 5 — IPS (Immunophenogram Score)")
            from row5_ips import IPSScorer
            IPSScorer().run(expr=expr_tpm, output_dir=str(self.output_dir))
        else:
            self.logger.info("Row 5 skipped (--skip-row5).")

        # ── Row 6: TIMER ───────────────────────────────────────────────────
        if not args.skip_row6:
            self.logger.info("=" * 60)
            self.logger.info("Row 6 — TIMER")
            if not args.cancer_type and not args.metadata:
                self.logger.error(
                    "Row 6 requires --cancer-type or --metadata. "
                    "Use --skip-row6 to bypass."
                )
                sys.exit(1)
            from row6_timer import TIMERScorer
            TIMERScorer().run(
                expr=expr_tpm,
                output_dir=str(self.output_dir),
                cancer_type=args.cancer_type,
                metadata_path=args.metadata,
            )
        else:
            self.logger.info("Row 6 skipped (--skip-row6).")

        # ── Row 7: Merge master table ──────────────────────────────────────
        self.logger.info("=" * 60)
        self.logger.info("Row 7 — Merging master table")
        from row7_merge import MasterMerger
        MasterMerger().run(
            output_dir=str(self.output_dir),
            timer_cancer_type=args.cancer_type,
        )

        # ── Row 8: Visualization ───────────────────────────────────────────
        if not args.skip_plots:
            self.logger.info("=" * 60)
            self.logger.info("Row 8 — Visualization")
            self._run_plots()
        else:
            self.logger.info("Row 8 skipped (--skip-plots).")

        self.logger.info("=" * 60)
        self.logger.info("Domain 2 pipeline complete. Outputs: %s", self.output_dir)

    # ── plot helper ───────────────────────────────────────────────────────

    def _run_plots(self) -> None:
        try:
            from plots import VisualizationSuite
        except ImportError as exc:
            self.logger.warning(
                "matplotlib/seaborn not available — skipping plots (%s). "
                "Install with: pip install matplotlib seaborn scikit-learn", exc,
            )
            return

        VisualizationSuite().run(
            output_dir=str(self.output_dir),
            cohort=self.args.cancer_type or "",
            dpi=self.args.dpi,
        )

    # ── clinical helper ───────────────────────────────────────────────────

    def _run_clinical(self) -> None:
        args = self.args

        if not args.clinical_patient or not args.clinical_sample:
            self.logger.info(
                "Row 0 — clinical harmonization skipped "
                "(pass --clinical-patient and --clinical-sample to enable)."
            )
            return

        self.logger.info("=" * 60)
        self.logger.info("Row 0 — Clinical Harmonization (%s)", args.clinical_source)

        from clinical_harmonizer import HARMONIZERS
        harmonizer_cls = HARMONIZERS.get(args.clinical_source)
        if harmonizer_cls is None:
            self.logger.error("Unknown --clinical-source '%s'.", args.clinical_source)
            return

        output_path = str(self.output_dir / "clinical_harmonized.csv")
        harmonizer_cls().run(
            output_path=output_path,
            patient_file=args.clinical_patient,
            sample_file=args.clinical_sample,
        )


def main() -> None:
    args = parse_args()
    _setup_logging(args.verbose)
    Domain2Pipeline(args).run()


if __name__ == "__main__":
    main()
