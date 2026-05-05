# Domain 2 — Immune Deconvolution & Microenvironment Profiling

Python/rpy2 pipeline for immune deconvolution of bulk RNA-seq expression data.  
Implements four complementary tools from the Domain 2 characterization layer:

| Row | Tool | Output | Value range |
|-----|------|--------|-------------|
| 3 | **EPIC · quanTIseq · MCP-counter** | Cell-type fraction matrix | 0–1 (fractions) / arbitrary units |
| 4 | **xCell** | 64 cell-type enrichment scores + 3 composite TME scores | 0–1 (enrichment) |
| 5 | **IPS** (Immunophenogram Score) | Score 0–10 + immunophenotype class A–D | 0–10 |
| 6 | **TIMER / TIMER2.0** | 6 TCGA-calibrated immune cell fractions | 0–1 (fractions) |

All outputs are merged into a single `Domain2_master_scores.csv` (samples × ~40 features).

---

## Requirements

### System
```bash
sudo apt-get install -y libmagick++-dev   # required by the GSVA → SpatialExperiment chain
```

### R (≥ 4.1)
```bash
export GITHUB_PAT=ghp_your_token_here    # needed — immunedeconv and several deps are GitHub-only
Rscript setup_r.R
```

`setup_r.R` installs all R dependencies in the correct order:
1. CRAN: `glmnet`, `igraph`, `magick`, `ComICS`, `DiagrammeR`, and immunedeconv core deps
2. Bioconductor: `BiocParallel`, `preprocessCore`, `Biobase`, `biomaRt`, `sva`, `quantiseqr`, `SpatialExperiment`, `GSVA`
3. GitHub-only: `EPIC`, `xCell`, `MCPcounter`, `mMCPcounter`, `ConsensusTME`
4. GitHub: `omnideconv/immunedeconv` (never use `BiocManager::install` for this)

### Python (≥ 3.10)
```bash
pip install -r requirements.txt
```

---

## Input Format

The pipeline accepts expression matrices in any of these forms — all handled automatically:

| Scenario | Handled automatically |
|----------|-----------------------|
| `genes × samples` (standard) | ✓ used as-is |
| `samples × genes` (transposed) | ✓ detected and flipped |
| Ensembl IDs (`ENSG...`) as gene names | ✓ mapped to HGNC symbols via MyGene.info |
| log2-transformed values (e.g. `log2(TPM+0.001)`) | ✓ offset parsed from filename and back-transformed |
| `.csv` or `.tsv` / `.txt` | ✓ separator auto-detected from extension |

**Required after loading:** linear TPM values, HGNC gene symbols as row index, sample IDs as column names.

---

## Usage

```bash
python main.py \
  --input  expression_tpm.tsv \
  --output ./results \
  --cancer-type SKCM
```

### All options

```
--input FILE           Expression matrix (CSV or TSV)
--output DIR           Output directory (default: ./results)
--cancer-type TYPE     TCGA cancer type for TIMER (e.g. SKCM, BRCA, LUAD)
--metadata CSV         Sample metadata with 'cancer_type' column (multi-cancer TIMER)
--log2-offset FLOAT    Pseudocount used before log2 (e.g. 0.001); auto-parsed from filename
--methods METHOD ...   Row 3 methods to run (default: epic quantiseq mcp_counter)
--cibersort-binary PATH  Path to CIBERSORT.R (optional, enables CIBERSORTx)
--cibersort-mat PATH     Path to LM22.txt (required with --cibersort-binary)
--skip-row3            Skip EPIC / quanTIseq / MCP-counter
--skip-row4            Skip xCell
--skip-row5            Skip IPS
--skip-row6            Skip TIMER
--verbose              Enable DEBUG-level logging
```

### Example — multi-cancer cohort

```bash
python main.py \
  --input  TCGA_HUGO_NORM_log2_TPM+0.001_mapped_matrix.tsv \
  --output ./results \
  --metadata metadata.csv    # must have columns: sample_id, cancer_type
```

---

## Output Files

All files are written to `--output` directory.

| File | Dimensions | Description |
|------|-----------|-------------|
| `deconv_epic.csv` | 8 × n | EPIC cell-type fractions (0–1) |
| `deconv_quantiseq.csv` | 11 × n | quanTIseq fractions (0–1, sum ≈ 1) |
| `deconv_mcp_counter.csv` | 10 × n | MCP-counter arbitrary units |
| `xcell_cell_scores.csv` | ≤64 × n | xCell per-cell-type enrichment scores |
| `xcell_composite_scores.csv` | 3 × n | ImmuneScore / StromaScore / MicroenvironmentScore |
| `ips_scores.csv` | n × 6 | IPS total, 4 module scores (MHC/EC/CP/SC), class A–D |
| `timer_{TYPE}.csv` | 6 × n | TIMER fractions calibrated to TCGA cancer type |
| `Domain2_master_scores.csv` | n × ~40 | All scores merged, samples as rows |

### Master table column prefixes

| Prefix | Source | Value type |
|--------|--------|------------|
| `EPIC_` | EPIC | Fraction 0–1 |
| `qT_` | quanTIseq | Fraction 0–1 |
| `MCP_` | MCP-counter | Arbitrary units |
| `xCell_` | xCell composite | Enrichment 0–1 |
| `IPS` / `IPS_class` | IPS | Score 0–10 / class A–D |
| `TIMER_` | TIMER | TCGA-calibrated fraction |

---

## Module Reference

| File | Class | Role |
|------|-------|------|
| `rpy2_bridge.py` | — | Shared pandas ↔ R matrix conversion; activates rpy2 converters |
| `input_validation.py` | — | Load, orient, map IDs, back-transform, validate |
| `row3_immune_fractions.py` | `ImmuneFrections` | EPIC / quanTIseq / MCP-counter / CIBERSORTx |
| `row4_xcell.py` | `XCellScorer` | xCell 64 cell types + composite scores |
| `row5_ips.py` | `IPSScorer` | IPS scoring + A–D classification |
| `row6_timer.py` | `TIMERScorer` | TIMER single- and multi-cancer-type scoring |
| `row7_merge.py` | `MasterMerger` | Merge all CSVs → master table + validation |
| `main.py` | `Domain2Pipeline` | CLI orchestrator |

---

## Validation Checks (automatic)

Run automatically after merging the master table:

- No NaN values anywhere in the master table
- EPIC fractions ≤ 1.05 per sample
- IPS scores within 0–10 for all samples
- TIMER fractions non-negative
- xCell `MicroenvironmentScore` ≥ 0
- EPIC vs quanTIseq CD8+ T cell Spearman ρ > 0.6 (logged as warning if not met)

---

## Important Notes

- **immunedeconv is not on Bioconductor** — `BiocManager::install("immunedeconv")` will fail. Use `setup_r.R` which installs from GitHub.
- **IPS requires log2 input** — `row5_ips.py` applies `log2(TPM+1)` internally. Always pass linear TPM to the pipeline.
- **TIMER is cancer-type specific** — using the wrong cancer type produces invalid results. Do not mix TIMER2.0 web portal scores with results from this pipeline.
- **xCell scores are not fractions** — they do not sum to 1 and cannot be compared to EPIC/quanTIseq values.
- **MCP-counter values are arbitrary units** — do not compare them to fraction-based methods.
