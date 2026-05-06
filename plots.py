"""
plots.py — Domain 2 immune deconvolution visualization suite.

Can be run standalone:
    python plots.py --results ./results [--dpi 150] [--cohort BRCA]

Or imported and called from the pipeline:
    from plots import VisualizationSuite
    VisualizationSuite().run(output_dir="./results", cohort="BRCA")

Figures produced in {output_dir}/plots/:
  1. epic_composition.png       EPIC cell fractions stacked bar (sorted by CD8+)
  2. quantiseq_composition.png  quanTIseq fractions stacked bar (same sample order)
  3. ips_profile.png            IPS score distribution + class bar chart
  4. xcell_composite.png        xCell composite TME scores by IPS class
  5. timer_fractions.png        TIMER 6-cell type violin plots
  6. cross_method_cd8.png       EPIC vs quanTIseq vs TIMER CD8+ scatter grid
  7. master_heatmap.png         Ward-clustered z-score heatmap of all numeric features

If clinical_harmonized.csv is present in output_dir, figures 3 and 7 gain
additional annotation tracks (stage, vital status).
"""
import argparse
import logging
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)
sns.set_theme(style="whitegrid", font_scale=1.05)

logger = logging.getLogger(__name__)

PALETTE_CLASS = {"A": "#2196F3", "B": "#4CAF50", "C": "#FF9800", "D": "#F44336"}
PALETTE_STAGE = {
    "I": "#E3F2FD", "IA": "#BBDEFB", "IB": "#90CAF9",
    "II": "#1565C0", "IIA": "#1976D2", "IIB": "#1E88E5",
    "III": "#880E4F", "IIIA": "#AD1457", "IIIB": "#C2185B", "IIIC": "#E91E63",
    "IV": "#B71C1C",
}


# ── helpers ───────────────────────────────────────────────────────────────────

def _save(fig, path: Path, dpi: int = 150) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("  saved → %s", path)


def _cohort_label(master: pd.DataFrame, cohort: str) -> str:
    """Derive cohort label from data if not supplied."""
    if cohort:
        return cohort
    if "cancer_type" in master.columns:
        ct = master["cancer_type"].dropna()
        if len(ct):
            return str(ct.iloc[0])
    return f"n={len(master)}"


def _stacked_bar(df_frac: pd.DataFrame, title: str, palette) -> plt.Figure:
    n = len(df_frac)
    fig, ax = plt.subplots(figsize=(max(14, n * 0.012), 5))
    bottom = np.zeros(n)
    for col, color in zip(df_frac.columns, palette):
        ax.bar(np.arange(n), df_frac[col].values, bottom=bottom,
               color=color, width=1.0, linewidth=0)
        bottom += df_frac[col].values
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel(f"Samples (n={n}, sorted by CD8+ T cell fraction)", labelpad=6)
    ax.set_ylabel("Fraction")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xticks([])
    handles = [mpatches.Patch(color=c, label=l)
               for l, c in zip(df_frac.columns, palette)]
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
              frameon=False, fontsize=9)
    return fig


# ── Figure 1 & 2: stacked bar compositions ───────────────────────────────────

def plot_compositions(master: pd.DataFrame, out: Path, cohort: str = "") -> None:
    label = _cohort_label(master, cohort)

    epic_cols = [c for c in master.columns if c.startswith("EPIC_")]
    if not epic_cols:
        logger.warning("No EPIC columns found — skipping composition plots.")
        return

    epic = master[epic_cols].rename(columns=lambda c: c.replace("EPIC_", ""))
    cd8_col = next((c for c in epic.columns if "CD8" in c), epic.columns[0])
    epic_sorted = epic.sort_values(cd8_col, ascending=False)

    fig = _stacked_bar(
        epic_sorted,
        f"EPIC  —  Cell-type fraction per sample  ({label})",
        sns.color_palette("tab10", len(epic.columns)),
    )
    _save(fig, out / "epic_composition.png")

    qt_cols = [c for c in master.columns if c.startswith("qT_")]
    if qt_cols:
        qt = master[qt_cols].rename(columns=lambda c: c.replace("qT_", ""))
        qt_sorted = qt.loc[epic_sorted.index]
        fig = _stacked_bar(
            qt_sorted,
            f"quanTIseq  —  Cell-type fraction per sample  (same order as EPIC)",
            sns.color_palette("tab20", len(qt.columns)),
        )
        _save(fig, out / "quantiseq_composition.png")


# ── Figure 3: IPS profile ─────────────────────────────────────────────────────

def plot_ips(master: pd.DataFrame, out: Path, cohort: str = "") -> None:
    if "IPS" not in master.columns or "IPS_class" not in master.columns:
        logger.warning("IPS columns not found — skipping IPS plot.")
        return

    label = _cohort_label(master, cohort)
    has_stage = "stage_overall" in master.columns and master["stage_overall"].notna().any()

    ncols = 3 if has_stage else 2
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
    fig.suptitle(f"Immunophenogram Score (IPS)  —  {label}", fontsize=13, fontweight="bold")

    # left: IPS histogram by class
    ax = axes[0]
    for cls, grp in master.groupby("IPS_class"):
        ax.hist(grp["IPS"], bins=20, color=PALETTE_CLASS.get(str(cls), "#999"),
                alpha=0.75, edgecolor="white", linewidth=0.4,
                label=f"Class {cls}  (n={len(grp)})")
    ax.set_xlabel("IPS score (0–10)")
    ax.set_ylabel("Sample count")
    ax.set_title("IPS distribution by class")
    ax.legend(frameon=False)

    # middle: class bar chart
    ax2 = axes[1]
    counts = master["IPS_class"].value_counts().sort_index()
    bars = ax2.bar(counts.index.astype(str), counts.values,
                   color=[PALETTE_CLASS.get(str(c), "#999") for c in counts.index],
                   edgecolor="white", linewidth=0.8)
    for bar, val in zip(bars, counts.values):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                 str(val), ha="center", va="bottom", fontsize=10)
    ax2.set_xlabel("IPS class")
    ax2.set_ylabel("Number of samples")
    ax2.set_title("Class distribution  (A≥7, B≥5, C≥3, D<3)")
    ax2.set_ylim(0, counts.max() * 1.15)

    # right (clinical): IPS by stage (only if clinical data available)
    if has_stage:
        ax3 = axes[2]
        stage_order = ["I", "IA", "IB", "II", "IIA", "IIB",
                       "III", "IIIA", "IIIB", "IIIC", "IV"]
        present = [s for s in stage_order if s in master["stage_overall"].values]
        pal = [PALETTE_STAGE.get(s, "#999") for s in present]
        sns.boxplot(data=master, x="stage_overall", y="IPS", order=present,
                    palette=pal, ax=ax3, linewidth=0.8, fliersize=2)
        ax3.set_xlabel("AJCC stage")
        ax3.set_ylabel("IPS score")
        ax3.set_title("IPS score by pathological stage")
        ax3.tick_params(axis="x", rotation=30)

    plt.tight_layout()
    _save(fig, out / "ips_profile.png")


# ── Figure 4: xCell composite scores ─────────────────────────────────────────

def plot_xcell(master: pd.DataFrame, out: Path, cohort: str = "") -> None:
    xcell_cols = [c for c in master.columns if c.startswith("xCell_")]
    if not xcell_cols or "IPS_class" not in master.columns:
        logger.warning("xCell or IPS_class columns not found — skipping xCell plot.")
        return

    df = master[xcell_cols + ["IPS_class"]].copy()
    score_cols = [c.replace("xCell_", "").title() for c in xcell_cols]
    df.columns = score_cols + ["IPS_class"]
    df_long = df.melt(id_vars="IPS_class", value_vars=score_cols,
                      var_name="Score", value_name="Value")

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.violinplot(data=df_long, x="Score", y="Value", hue="IPS_class",
                   palette=PALETTE_CLASS, inner="box", cut=0, ax=ax,
                   hue_order=sorted(master["IPS_class"].unique()))
    ax.set_title("xCell Composite TME Scores  —  by IPS class",
                 fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("xCell enrichment score")
    ax.legend(title="IPS class", bbox_to_anchor=(1.01, 1), frameon=False)
    plt.tight_layout()
    _save(fig, out / "xcell_composite.png")


# ── Figure 5: TIMER cell fractions ────────────────────────────────────────────

def plot_timer(master: pd.DataFrame, out: Path, cohort: str = "") -> None:
    timer_cols = [c for c in master.columns if c.startswith("TIMER_")]
    if not timer_cols:
        logger.warning("No TIMER columns found — skipping TIMER plot.")
        return

    label = _cohort_label(master, cohort)
    df = master[timer_cols].rename(columns=lambda c: c.replace("TIMER_", ""))
    order = df.median().sort_values(ascending=False).index.tolist()

    fig, ax = plt.subplots(figsize=(11, 5))
    sns.violinplot(data=df.melt(var_name="Cell type", value_name="Fraction"),
                   x="Cell type", y="Fraction", order=order,
                   palette="Blues_d", inner="box", cut=0, ax=ax)
    ax.set_title(f"TIMER  —  Immune cell fractions  ({label}, TCGA-calibrated)",
                 fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Estimated fraction")
    ax.tick_params(axis="x", rotation=20)
    plt.tight_layout()
    _save(fig, out / "timer_fractions.png")


# ── Figure 6: cross-method CD8+ scatter ──────────────────────────────────────

def plot_cross_method_cd8(master: pd.DataFrame, out: Path) -> None:
    needed = {"EPIC_T cell CD8+", "qT_T cell CD8+", "TIMER_T cell CD8+"}
    if not needed.issubset(master.columns):
        logger.warning("Missing CD8+ columns for cross-method plot — skipping.")
        return

    epic  = master["EPIC_T cell CD8+"]
    qt    = master["qT_T cell CD8+"]
    timer = master["TIMER_T cell CD8+"]

    pairs = [
        (epic, qt,    "EPIC",       "quanTIseq"),
        (epic, timer, "EPIC",       "TIMER"),
        (qt,   timer, "quanTIseq",  "TIMER"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Cross-method CD8+ T cell fraction comparison",
                 fontsize=13, fontweight="bold")

    for ax, (x, y, xlabel, ylabel) in zip(axes, pairs):
        mask = (x > 0) | (y > 0)
        rho, pval = stats.spearmanr(x[mask], y[mask])
        ax.scatter(x[mask], y[mask], s=5, alpha=0.35, color="#1565C0", rasterized=True)
        lim = max(x[mask].max(), y[mask].max()) * 1.05
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.4)
        ax.set_xlim(-0.005, lim)
        ax.set_ylim(-0.005, lim)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(f"ρ = {rho:.2f}  (p={pval:.1e}, n={mask.sum()})", fontsize=10)

    plt.tight_layout()
    _save(fig, out / "cross_method_cd8.png")


# ── Figure 7: master feature heatmap ─────────────────────────────────────────

def plot_master_heatmap(master: pd.DataFrame, out: Path, cohort: str = "",
                        dpi: int = 120) -> None:
    numeric = master.select_dtypes(include="number").copy()
    # Drop columns with >20% missing, then median-impute any residual NaN
    numeric = numeric.loc[:, numeric.isna().mean() < 0.2]
    numeric = numeric.fillna(numeric.median())
    if numeric.shape[1] < 2 or not np.isfinite(numeric.values).all():
        logger.warning("Not enough finite numeric columns for heatmap — skipping.")
        return

    label = _cohort_label(master, cohort)

    # z-score, clip ±3
    z = pd.DataFrame(
        StandardScaler().fit_transform(numeric),
        index=numeric.index, columns=numeric.columns,
    ).clip(-3, 3)

    var_order = z.var().sort_values(ascending=False).index[:min(300, len(z.columns))]
    z_sub = z[var_order]

    row_order = dendrogram(linkage(z_sub.values,   method="ward"), no_plot=True)["leaves"]
    col_order = dendrogram(linkage(z_sub.values.T, method="ward"), no_plot=True)["leaves"]
    z_plot = z_sub.iloc[row_order, col_order]

    # Determine how many annotation rows to show
    ann_tracks: list[tuple[str, np.ndarray]] = []

    # Track 1: IPS class (always available if IPS ran)
    if "IPS_class" in master.columns:
        ips_classes = master.loc[z_plot.index, "IPS_class"].astype(str)
        rgba = np.array([mcolors.to_rgba(PALETTE_CLASS.get(c, "#999999"))
                         for c in ips_classes], dtype=float)
        ann_tracks.append(("IPS class", rgba))

    # Track 2: stage (clinical data)
    if "stage_overall" in master.columns:
        stages = master.loc[z_plot.index, "stage_overall"].fillna("unknown").astype(str)
        all_stages = sorted(stages.unique())
        stage_cmap = plt.cm.get_cmap("tab20", len(all_stages))
        stage_idx = {s: i for i, s in enumerate(all_stages)}
        rgba_stage = np.array([mcolors.to_rgba(stage_cmap(stage_idx[s]))
                                for s in stages], dtype=float)
        ann_tracks.append(("Stage", rgba_stage))

    # Track 3: vital status (clinical data)
    if "vital_status" in master.columns:
        vs = master.loc[z_plot.index, "vital_status"].fillna("unknown").astype(str)
        vs_palette = {"alive": "#43A047", "dead": "#E53935", "unknown": "#BDBDBD"}
        rgba_vs = np.array([mcolors.to_rgba(vs_palette.get(v, "#BDBDBD"))
                            for v in vs], dtype=float)
        ann_tracks.append(("Vital status", rgba_vs))

    n_ann = len(ann_tracks)
    height_ratios = [0.022] * n_ann + [1]
    fig, axes = plt.subplots(
        n_ann + 1, 1,
        figsize=(18, 12),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.01},
    )
    if n_ann == 0:
        axes = [axes]

    for i, (track_label, rgba) in enumerate(ann_tracks):
        ax = axes[i]
        ax.imshow(rgba.reshape(1, len(rgba), 4), aspect="auto", interpolation="none")
        ax.set_xticks([])
        ax.set_yticks([0])
        ax.set_yticklabels([track_label], fontsize=8)

    hm_ax = axes[n_ann]
    im = hm_ax.imshow(z_plot.values.T, aspect="auto", cmap="RdBu_r",
                      vmin=-3, vmax=3, interpolation="none")
    hm_ax.set_xticks([])
    hm_ax.set_yticks(range(len(z_plot.columns)))
    hm_ax.set_yticklabels(z_plot.columns, fontsize=7)
    hm_ax.set_xlabel(f"Samples (n={len(z_plot)}, Ward-clustered)", fontsize=10)
    hm_ax.set_title(
        f"Domain 2 master feature heatmap  —  z-scored, Ward-clustered  ({label})",
        fontsize=13, fontweight="bold", pad=4,
    )
    fig.colorbar(im, ax=hm_ax, orientation="vertical",
                 fraction=0.015, pad=0.01).set_label("z-score (clipped ±3)", fontsize=8)

    # Legends
    legend_handles = []
    if "IPS_class" in master.columns:
        legend_handles += [mpatches.Patch(color=v, label=f"IPS {k}")
                           for k, v in PALETTE_CLASS.items()
                           if k in master["IPS_class"].unique().astype(str)]
    if "vital_status" in master.columns:
        for label_vs, col_vs in {"alive": "#43A047", "dead": "#E53935"}.items():
            legend_handles.append(mpatches.Patch(color=col_vs, label=label_vs))
    if legend_handles:
        fig.legend(handles=legend_handles, loc="lower right",
                   frameon=False, fontsize=9, ncol=2)

    _save(fig, out / "master_heatmap.png", dpi=dpi)


# ── public entry point (used by pipeline) ─────────────────────────────────────

class VisualizationSuite:
    """Run all Domain 2 visualizations from a results directory."""

    def run(self, output_dir: str, cohort: str = "", dpi: int = 150) -> Path:
        """
        Parameters
        ----------
        output_dir : directory containing Domain2_master_scores.csv
                     (and optionally clinical_harmonized.csv)
        cohort     : label for plot titles (e.g. 'BRCA'); auto-derived if empty
        dpi        : figure resolution

        Returns
        -------
        Path to the plots sub-directory
        """
        out_dir = Path(output_dir)
        plots_dir = out_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        # Prefer master+clinical when available
        master_clinical = out_dir / "Domain2_master_with_clinical.csv"
        master_only     = out_dir / "Domain2_master_scores.csv"

        if master_clinical.exists():
            master = pd.read_csv(master_clinical, index_col=0, low_memory=False)
            logger.info("Plotting from master+clinical table (%d samples × %d features)",
                        master.shape[0], master.shape[1])
        elif master_only.exists():
            master = pd.read_csv(master_only, index_col=0)
            logger.info("Plotting from deconvolution master table (%d samples × %d features)",
                        master.shape[0], master.shape[1])
        else:
            logger.error("No master table found in %s — run the pipeline first.", out_dir)
            return plots_dir

        steps = [
            ("EPIC & quanTIseq composition",  lambda: plot_compositions(master, plots_dir, cohort)),
            ("IPS profile",                   lambda: plot_ips(master, plots_dir, cohort)),
            ("xCell composite scores",        lambda: plot_xcell(master, plots_dir, cohort)),
            ("TIMER fractions",               lambda: plot_timer(master, plots_dir, cohort)),
            ("Cross-method CD8+ comparison",  lambda: plot_cross_method_cd8(master, plots_dir)),
            ("Master feature heatmap",        lambda: plot_master_heatmap(master, plots_dir, cohort, dpi)),
        ]

        for i, (name, fn) in enumerate(steps, 1):
            logger.info("[%d/%d] %s…", i, len(steps), name)
            try:
                fn()
            except Exception as exc:
                logger.warning("  %s failed: %s", name, exc)

        logger.info("Plots complete → %s", plots_dir)
        return plots_dir


# ── standalone CLI ─────────────────────────────────────────────────────────────

def _cli() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)-8s %(message)s")
    ap = argparse.ArgumentParser(description="Domain 2 visualization suite")
    ap.add_argument("--results", default="results",
                    help="Results directory (default: results)")
    ap.add_argument("--cohort",  default="",
                    help="Cohort label for plot titles (auto-detected if omitted)")
    ap.add_argument("--dpi",    type=int, default=150)
    args = ap.parse_args()
    VisualizationSuite().run(output_dir=args.results, cohort=args.cohort, dpi=args.dpi)


if __name__ == "__main__":
    _cli()
