"""
input_validation.py — load and validate the expression matrix before deconvolution.

Handles three common real-world issues automatically:
  1. Transposed input (samples×genes) — detected and flipped to genes×samples
  2. Ensembl IDs (ENSG...) — mapped to HGNC symbols via MyGene.info
  3. log2-transformed values — back-transformed to linear TPM using the
     pseudocount parsed from the filename (e.g. "log2_TPM+0.001" → 0.001)
     or supplied explicitly via log2_offset.
"""
import logging
import re
from pathlib import Path

import mygene
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

_LOG2_OFFSET_RE = re.compile(r"log2[_\-]?tpm\+?([\d.]+)", re.IGNORECASE)


# ── public API ────────────────────────────────────────────────────────────────

def load_expression_matrix(
    path: str,
    log2_offset: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load a genes×samples expression matrix from a CSV or TSV file.

    Parameters
    ----------
    path        : path to the file (.tsv/.txt → tab-separated, else comma)
    log2_offset : pseudocount that was added before log2 transformation
                  (e.g. 0.001 for log2(TPM+0.001)).  Pass None to skip
                  back-transformation; the function also auto-detects the
                  offset from the filename when log2 is present.

    Returns
    -------
    expr_tpm  : linear TPM DataFrame, genes×samples, HGNC row index
    expr_log2 : log2(expr_tpm + 1) — used internally by IPS
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Expression matrix not found: {path}")

    # ── 1. Read file ──────────────────────────────────────────────────────────
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    expr = pd.read_csv(path, index_col=0, sep=sep)

    if expr.shape[0] == 0 or expr.shape[1] == 0:
        raise ValueError(
            f"Matrix is empty after reading '{path.name}' with sep={repr(sep)}."
        )

    # ── 2. Orient to genes×samples ────────────────────────────────────────────
    expr = _ensure_genes_as_rows(expr)

    # ── 3. Map Ensembl IDs → HGNC symbols ────────────────────────────────────
    if _looks_like_ensembl(expr.index):
        logger.info("Ensembl IDs detected — mapping to HGNC symbols via MyGene.info...")
        expr = _map_ensembl_to_hgnc(expr)

    # ── 4. Back-transform log2 → linear TPM ──────────────────────────────────
    offset = log2_offset if log2_offset is not None else _parse_log2_offset(path.name)
    if offset is not None:
        logger.info(
            "log2 input detected (offset=%.4g) — back-transforming to linear TPM: "
            "2^x − %.4g, values clipped at 0.",
            offset, offset,
        )
        idx, cols = expr.index, expr.columns
        expr = pd.DataFrame(
            np.maximum(np.power(2.0, expr.values) - offset, 0.0),
            index=idx,
            columns=cols,
        )
    elif float(expr.values.max()) < 30:
        logger.warning(
            "Max expression value is %.2f — matrix may be log-transformed. "
            "If so, pass --log2-offset (e.g. --log2-offset 1) to back-transform. "
            "EPIC/quanTIseq/xCell/TIMER require linear TPM.",
            float(expr.values.max()),
        )

    # ── 5. Final validation ───────────────────────────────────────────────────
    if expr.isnull().sum().sum() != 0:
        raise ValueError(f"Matrix contains {expr.isnull().sum().sum()} NaN values.")
    if not (expr >= 0).all().all():
        raise ValueError("Matrix contains negative values — TPM values must be >= 0.")

    logger.info(
        "Expression matrix ready: %d genes × %d samples  "
        "(TPM range: %.2f – %.2f)",
        expr.shape[0], expr.shape[1],
        float(expr.values.min()), float(expr.values.max()),
    )

    expr_log2 = np.log2(expr + 1)
    return expr, expr_log2


# ── helpers ───────────────────────────────────────────────────────────────────

def _looks_like_ensembl(index: pd.Index) -> bool:
    """Return True if >50% of index entries look like Ensembl gene IDs."""
    sample = index[:min(20, len(index))]
    hits = sum(1 for s in sample if str(s).startswith("ENSG"))
    return hits / len(sample) > 0.5


def _looks_like_samples(index: pd.Index) -> bool:
    """Return True if index entries look like sample IDs (not gene symbols)."""
    sample = index[:min(20, len(index))]
    # TCGA-style: "TCGA-AB-1234-01", or plain alphanumeric with hyphens > 2 parts
    hits = sum(
        1 for s in sample
        if str(s).startswith("TCGA-") or len(str(s).split("-")) >= 3
    )
    return hits / len(sample) > 0.5


def _ensure_genes_as_rows(expr: pd.DataFrame) -> pd.DataFrame:
    """Transpose if rows look like sample IDs rather than gene symbols."""
    if _looks_like_samples(expr.index) or _looks_like_ensembl(expr.columns):
        logger.info(
            "Input appears to be samples×genes (%d × %d) — transposing to genes×samples.",
            expr.shape[0], expr.shape[1],
        )
        expr = expr.T
    return expr


def _map_ensembl_to_hgnc(expr: pd.DataFrame) -> pd.DataFrame:
    """Map Ensembl row index to HGNC gene symbols; drop unmapped genes."""
    mg_client = mygene.MyGeneInfo()
    ensembl_ids = expr.index.tolist()

    results = mg_client.querymany(
        ensembl_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=True,
        df_index=True,
    )

    # Build Ensembl → symbol mapping (drop duplicates and missing)
    sym_col = "symbol"
    if sym_col not in results.columns:
        raise RuntimeError("MyGene.info returned no 'symbol' column — check network access.")

    mapping = results[[sym_col]].dropna()
    mapping = mapping[~mapping.index.duplicated(keep=False)]      # drop ambiguous Ensembl IDs
    mapping = mapping[~mapping[sym_col].duplicated(keep=False)]   # drop ambiguous symbols

    # Remap index — both index and symbol are now 1-to-1, so .loc is safe
    expr_mapped = expr.loc[expr.index.isin(mapping.index)].copy()
    expr_mapped.index = mapping.loc[expr_mapped.index, sym_col].values

    # Drop any duplicate HGNC symbols (keep mean)
    if expr_mapped.index.duplicated().any():
        n_dup = expr_mapped.index.duplicated().sum()
        logger.warning("Collapsing %d duplicate HGNC symbols by mean.", n_dup)
        expr_mapped = expr_mapped.groupby(expr_mapped.index).mean()

    n_mapped = len(expr_mapped)
    n_total = len(expr)
    logger.info("Ensembl → HGNC: %d / %d genes mapped.", n_mapped, n_total)
    if n_mapped < n_total * 0.5:
        logger.warning(
            "Fewer than 50%% of genes mapped (%d / %d). "
            "Verify the gene IDs are human Ensembl gene IDs (ENSG...).",
            n_mapped, n_total,
        )
    return expr_mapped


def _parse_log2_offset(filename: str) -> float | None:
    """
    Extract the pseudocount from filenames like 'log2_TPM+0.001' or 'log2TPM+1'.
    Returns None if no log2 pattern is found.
    """
    m = _LOG2_OFFSET_RE.search(filename)
    if m:
        return float(m.group(1))
    # filename contains "log2" but no explicit offset — assume +1
    if re.search(r"log2", filename, re.IGNORECASE):
        logger.info("'log2' in filename but no offset found — assuming pseudocount = 1.")
        return 1.0
    return None
