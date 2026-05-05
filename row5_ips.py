"""
row5_ips.py — Row 5: IPS (Immunophenogram Score)
Reference: Charoentong et al. 2017 Cell Reports 18:248-262
Gene table: icbi-lab/Immunophenogram (IPS_genes.txt)

immunedeconv 2.x (omnideconv) removed ips_rand entirely.
This module reimplements the original ICBI algorithm directly in Python,
matching the R reference code from icbi-lab/Immunophenogram/IPS.R exactly.

Algorithm (per sample):
  1. Z-score every IPS gene using the GLOBAL sample mean and SD (all genes).
  2. For each of the 26 unique gene groups (NAME column), average the z-scores
     of all genes in that group, then multiply by the group weight.
  3. Aggregate: MHC = mean(groups 1-10), CP = mean(groups 11-20),
                EC  = mean(groups 21-24), SC = mean(groups 25-26)
  4. AZ = MHC + CP + EC + SC
  5. IPS = 0 if AZ≤0, 10 if AZ≥3, else round(AZ*10/3)

Input: log2(TPM+1) — applied internally from linear TPM.
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ── IPS gene table (icbi-lab/Immunophenogram/IPS_genes.txt) ─────────────────
# Columns: GENE (HGNC symbol), NAME (group label), CLASS (MHC/EC/CP/SC), WEIGHT
_IPS_GENES = [
    # MHC — groups 1-10
    ("B2M",       "B2M",       "MHC",  1),
    ("TAP1",      "TAP1",      "MHC",  1),
    ("TAP2",      "TAP2",      "MHC",  1),
    ("HLA-A",     "HLA-A",     "MHC",  1),
    ("HLA-B",     "HLA-B",     "MHC",  1),
    ("HLA-C",     "HLA-C",     "MHC",  1),
    ("HLA-DPA1",  "HLA-DPA1",  "MHC",  1),
    ("HLA-DPB1",  "HLA-DPB1",  "MHC",  1),
    ("HLA-E",     "HLA-E",     "MHC",  1),
    ("HLA-F",     "HLA-F",     "MHC",  1),
    # CP — groups 11-20
    ("PDCD1",     "PD-1",      "CP",  -1),
    ("CTLA4",     "CTLA-4",    "CP",  -1),
    ("LAG3",      "LAG3",      "CP",  -1),
    ("TIGIT",     "TIGIT",     "CP",  -1),
    ("HAVCR2",    "TIM3",      "CP",  -1),
    ("CD274",     "PD-L1",     "CP",  -1),
    ("PDCD1LG2",  "PD-L2",     "CP",  -1),
    ("CD27",      "CD27",      "CP",   1),
    ("ICOS",      "ICOS",      "CP",   1),
    ("IDO1",      "IDO1",      "CP",  -1),
    # EC — groups 21-24 (multiple genes share same NAME → averaged)
    ("AIM2",      "Act CD4",   "EC",  1), ("BIRC3",   "Act CD4",  "EC",  1),
    ("BRIP1",     "Act CD4",   "EC",  1), ("CCL20",   "Act CD4",  "EC",  1),
    ("CCL4",      "Act CD4",   "EC",  1), ("CCL5",    "Act CD4",  "EC",  1),
    ("CCNB1",     "Act CD4",   "EC",  1), ("CCR7",    "Act CD4",  "EC",  1),
    ("DUSP2",     "Act CD4",   "EC",  1), ("ESCO2",   "Act CD4",  "EC",  1),
    ("ETS1",      "Act CD4",   "EC",  1), ("EXO1",    "Act CD4",  "EC",  1),
    ("EXOC6",     "Act CD4",   "EC",  1), ("IARS",    "Act CD4",  "EC",  1),
    ("KIF11",     "Act CD4",   "EC",  1), ("KNTC1",   "Act CD4",  "EC",  1),
    ("NUF2",      "Act CD4",   "EC",  1), ("PRC1",    "Act CD4",  "EC",  1),
    ("PSAT1",     "Act CD4",   "EC",  1), ("RGS1",    "Act CD4",  "EC",  1),
    ("RTKN2",     "Act CD4",   "EC",  1), ("SAMSN1",  "Act CD4",  "EC",  1),
    ("SELL",      "Act CD4",   "EC",  1), ("TRAT1",   "Act CD4",  "EC",  1),
    ("ADRM1",     "Act CD8",   "EC",  1), ("AHSA1",   "Act CD8",  "EC",  1),
    ("C1GALT1C1", "Act CD8",   "EC",  1), ("CCT6B",   "Act CD8",  "EC",  1),
    ("CD37",      "Act CD8",   "EC",  1), ("CD3D",    "Act CD8",  "EC",  1),
    ("CD3E",      "Act CD8",   "EC",  1), ("CD3G",    "Act CD8",  "EC",  1),
    ("CD69",      "Act CD8",   "EC",  1), ("CD8A",    "Act CD8",  "EC",  1),
    ("CETN3",     "Act CD8",   "EC",  1), ("CSE1L",   "Act CD8",  "EC",  1),
    ("GEMIN6",    "Act CD8",   "EC",  1), ("GNLY",    "Act CD8",  "EC",  1),
    ("GPT2",      "Act CD8",   "EC",  1), ("GZMA",    "Act CD8",  "EC",  1),
    ("GZMH",      "Act CD8",   "EC",  1), ("GZMK",    "Act CD8",  "EC",  1),
    ("IL2RB",     "Act CD8",   "EC",  1), ("LCK",     "Act CD8",  "EC",  1),
    ("MPZL1",     "Act CD8",   "EC",  1), ("NKG7",    "Act CD8",  "EC",  1),
    ("PIK3IP1",   "Act CD8",   "EC",  1), ("PTRH2",   "Act CD8",  "EC",  1),
    ("TIMM13",    "Act CD8",   "EC",  1), ("ZAP70",   "Act CD8",  "EC",  1),
    ("ATM",       "Tem CD4",   "EC",  1), ("CASP3",   "Tem CD4",  "EC",  1),
    ("CASQ1",     "Tem CD4",   "EC",  1), ("CD300E",  "Tem CD4",  "EC",  1),
    ("DARS",      "Tem CD4",   "EC",  1), ("DOCK9",   "Tem CD4",  "EC",  1),
    ("EXOSC9",    "Tem CD4",   "EC",  1), ("EZH2",    "Tem CD4",  "EC",  1),
    ("GDE1",      "Tem CD4",   "EC",  1), ("IL34",    "Tem CD4",  "EC",  1),
    ("NCOA4",     "Tem CD4",   "EC",  1), ("NEFL",    "Tem CD4",  "EC",  1),
    ("PDGFRL",    "Tem CD4",   "EC",  1), ("PTGS1",   "Tem CD4",  "EC",  1),
    ("REPS1",     "Tem CD4",   "EC",  1), ("SCG2",    "Tem CD4",  "EC",  1),
    ("SDPR",      "Tem CD4",   "EC",  1), ("SIGLEC14","Tem CD4",  "EC",  1),
    ("SIGLEC6",   "Tem CD4",   "EC",  1), ("TAL1",    "Tem CD4",  "EC",  1),
    ("TFEC",      "Tem CD4",   "EC",  1), ("TIPIN",   "Tem CD4",  "EC",  1),
    ("TPK1",      "Tem CD4",   "EC",  1), ("UQCRB",   "Tem CD4",  "EC",  1),
    ("USP9Y",     "Tem CD4",   "EC",  1), ("WIPF1",   "Tem CD4",  "EC",  1),
    ("ZCRB1",     "Tem CD4",   "EC",  1),
    ("ACAP1",     "Tem CD8",   "EC",  1), ("APOL3",   "Tem CD8",  "EC",  1),
    ("ARHGAP10",  "Tem CD8",   "EC",  1), ("ATP10D",  "Tem CD8",  "EC",  1),
    ("C3AR1",     "Tem CD8",   "EC",  1), ("CCR5",    "Tem CD8",  "EC",  1),
    ("CD160",     "Tem CD8",   "EC",  1), ("CD55",    "Tem CD8",  "EC",  1),
    ("CFLAR",     "Tem CD8",   "EC",  1), ("CMKLR1",  "Tem CD8",  "EC",  1),
    ("DAPP1",     "Tem CD8",   "EC",  1), ("FCRL6",   "Tem CD8",  "EC",  1),
    ("FLT3LG",    "Tem CD8",   "EC",  1), ("GZMM",    "Tem CD8",  "EC",  1),
    ("HAPLN3",    "Tem CD8",   "EC",  1), ("HLA-DMB", "Tem CD8",  "EC",  1),
    ("HLA-DPA1",  "Tem CD8",   "EC",  1), ("HLA-DPB1","Tem CD8",  "EC",  1),
    ("IFI16",     "Tem CD8",   "EC",  1), ("LIME1",   "Tem CD8",  "EC",  1),
    ("LTK",       "Tem CD8",   "EC",  1), ("NFKBIA",  "Tem CD8",  "EC",  1),
    ("SETD7",     "Tem CD8",   "EC",  1), ("SIK1",    "Tem CD8",  "EC",  1),
    ("TRIB2",     "Tem CD8",   "EC",  1),
    # SC — groups 25-26
    ("CCR2",      "MDSC",      "SC", -1), ("CD14",    "MDSC",     "SC", -1),
    ("CD2",       "MDSC",      "SC", -1), ("CD86",    "MDSC",     "SC", -1),
    ("CXCR4",     "MDSC",      "SC", -1), ("FCGR2A",  "MDSC",     "SC", -1),
    ("FCGR2B",    "MDSC",      "SC", -1), ("FCGR3A",  "MDSC",     "SC", -1),
    ("FERMT3",    "MDSC",      "SC", -1), ("GPSM3",   "MDSC",     "SC", -1),
    ("IL18BP",    "MDSC",      "SC", -1), ("IL4R",    "MDSC",     "SC", -1),
    ("ITGAL",     "MDSC",      "SC", -1), ("ITGAM",   "MDSC",     "SC", -1),
    ("PARVG",     "MDSC",      "SC", -1), ("PSAP",    "MDSC",     "SC", -1),
    ("PTGER2",    "MDSC",      "SC", -1), ("PTGES2",  "MDSC",     "SC", -1),
    ("S100A8",    "MDSC",      "SC", -1), ("S100A9",  "MDSC",     "SC", -1),
    ("CCL3L1",    "Treg",      "SC", -1), ("CD72",    "Treg",     "SC", -1),
    ("CLEC5A",    "Treg",      "SC", -1), ("FOXP3",   "Treg",     "SC", -1),
    ("ITGA4",     "Treg",      "SC", -1), ("L1CAM",   "Treg",     "SC", -1),
    ("LIPA",      "Treg",      "SC", -1), ("LRP1",    "Treg",     "SC", -1),
    ("LRRC42",    "Treg",      "SC", -1), ("MARCO",   "Treg",     "SC", -1),
    ("MMP12",     "Treg",      "SC", -1), ("MNDA",    "Treg",     "SC", -1),
    ("MRC1",      "Treg",      "SC", -1), ("MS4A6A",  "Treg",     "SC", -1),
    ("PELO",      "Treg",      "SC", -1), ("PLEK",    "Treg",     "SC", -1),
    ("PRSS23",    "Treg",      "SC", -1), ("PTGIR",   "Treg",     "SC", -1),
    ("ST8SIA4",   "Treg",      "SC", -1), ("STAB1",   "Treg",     "SC", -1),
]

_IPS_TABLE = pd.DataFrame(_IPS_GENES, columns=["GENE", "NAME", "CLASS", "WEIGHT"])

# Ordered list of the 26 unique group names (preserves file order, matches R code)
_UNIQUE_NAMES = list(dict.fromkeys(_IPS_TABLE["NAME"]))  # insertion-order dedup


def _ipsmap(az: float) -> int:
    """Map continuous AZ score to integer IPS 0-10 (matches R ipsmap function)."""
    if az <= 0:
        return 0
    if az >= 3:
        return 10
    return int(round(az * 10 / 3, 0))


def _classify(ips: int) -> str:
    if ips >= 7:
        return "A"
    if ips >= 5:
        return "B"
    if ips >= 3:
        return "C"
    return "D"


class IPSScorer:
    """Score samples with the Immunophenogram (IPS) and assign class A–D."""

    def run(self, expr: pd.DataFrame, output_dir: str) -> pd.DataFrame:
        """
        Parameters
        ----------
        expr       : LINEAR TPM DataFrame (genes × samples, HGNC row index).
                     log2(TPM+1) transform is applied here.
        output_dir : directory where ips_scores.csv is saved.

        Returns
        -------
        DataFrame indexed by sample ID with columns:
        IPS_score, MHC, CP, EC, SC, Class
        """
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        # IPS requires log2(TPM+1)
        expr_log2 = np.log2(expr + 1)
        logger.info("Running IPS on %d genes × %d samples (pure Python)...",
                    expr_log2.shape[0], expr_log2.shape[1])

        # Report gene coverage
        present = set(expr_log2.index) & set(_IPS_TABLE["GENE"])
        total = len(set(_IPS_TABLE["GENE"]))
        logger.info("IPS gene coverage: %d / %d genes present.", len(present), total)
        if len(present) < total * 0.6:
            logger.warning(
                "Fewer than 60%% of IPS genes found (%d / %d). "
                "Results may be unreliable.",
                len(present), total,
            )

        records = []
        for sample in expr_log2.columns:
            ge = expr_log2[sample]
            mge = ge.mean()
            sge = ge.std(ddof=1)

            # Z-scores for IPS genes using global sample mean/SD
            ips_genes_present = _IPS_TABLE["GENE"].isin(ge.index)
            tbl = _IPS_TABLE[ips_genes_present].copy()
            tbl = tbl[tbl["GENE"].isin(ge.index)].copy()
            tbl["Z"] = (ge.reindex(tbl["GENE"]).values - mge) / sge

            # Per-group weighted score: mean(z for genes in group) * group weight
            wg = {}
            for name in _UNIQUE_NAMES:
                mask = tbl["NAME"] == name
                if not mask.any():
                    wg[name] = 0.0
                    continue
                grp = tbl[mask]
                wg[name] = float(grp["Z"].mean()) * float(grp["WEIGHT"].iloc[0])

            mhc = float(np.mean([wg[n] for n in _UNIQUE_NAMES[0:10]]))
            cp  = float(np.mean([wg[n] for n in _UNIQUE_NAMES[10:20]]))
            ec  = float(np.mean([wg[n] for n in _UNIQUE_NAMES[20:24]]))
            sc  = float(np.mean([wg[n] for n in _UNIQUE_NAMES[24:26]]))
            az  = mhc + cp + ec + sc
            ips = _ipsmap(az)

            records.append({
                "IPS_score": ips,
                "MHC": round(mhc, 4),
                "CP":  round(cp,  4),
                "EC":  round(ec,  4),
                "SC":  round(sc,  4),
                "Class": _classify(ips),
            })

        result = pd.DataFrame(records, index=expr_log2.columns)
        result.to_csv(out / "ips_scores.csv")

        class_counts = result["Class"].value_counts().to_dict()
        logger.info("IPS done. %d samples. Classes: %s", len(result), class_counts)

        ips_min, ips_max = int(result["IPS_score"].min()), int(result["IPS_score"].max())
        if not (0 <= ips_min and ips_max <= 10):
            logger.warning("IPS scores outside 0-10 (min=%d, max=%d).", ips_min, ips_max)

        return result
