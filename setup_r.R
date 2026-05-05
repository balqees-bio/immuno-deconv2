# setup_r.R — install R dependencies for Domain 2 immune deconvolution
#
# immunedeconv is NOT on Bioconductor. Install via GitHub (omnideconv/immunedeconv).
# Several dependencies are also GitHub-only: EPIC, xCell, MCPcounter, ConsensusTME.
#
# BEFORE running this script:
#   1. Install system library (already done if you followed setup):
#        sudo apt-get install -y libmagick++-dev
#   2. Create a GitHub Personal Access Token (no special scopes needed):
#        https://github.com/settings/tokens  -> "Generate new token (classic)"
#   3. Export it then run:
#        export GITHUB_PAT=ghp_xxxxxxxxxxxx
#        Rscript setup_r.R

# ── 0. Require GITHUB_PAT — fail early with clear instructions ────────────────
pat <- Sys.getenv("GITHUB_PAT")
if (nchar(pat) == 0) {
  stop(
    "\nGITHUB_PAT is not set. GitHub throttles unauthenticated installs to\n",
    "60 requests/hour which is not enough for all dependencies.\n\n",
    "Fix:\n",
    "  1. Go to https://github.com/settings/tokens\n",
    "  2. Generate a new classic token (no scopes needed for public repos)\n",
    "  3. Run:  export GITHUB_PAT=ghp_your_token_here\n",
    "  4. Re-run: Rscript setup_r.R\n"
  )
}
cat("GITHUB_PAT is set.\n")

CRAN <- "https://cloud.r-project.org"
options(repos = c(CRAN = CRAN))

# ── 1. Remove stale lock directories from any previous failed run ─────────────
lib <- .libPaths()[1]
locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
if (length(locks) > 0) {
  cat("Removing stale locks:", paste(basename(locks), collapse = ", "), "\n")
  unlink(locks, recursive = TRUE)
}

# ── helpers ───────────────────────────────────────────────────────────────────
needs <- function(pkgs) pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

clear_locks <- function() {
  lk <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
  if (length(lk)) unlink(lk, recursive = TRUE)
}

gh_install <- function(repo) {
  pkg <- sub(".*/", "", repo)
  cat("Installing", pkg, "from GitHub (", repo, ")...\n")
  tryCatch(
    remotes::install_github(repo, upgrade = "never", force = TRUE),
    error = function(e) cat("  WARNING:", conditionMessage(e), "\n")
  )
  clear_locks()
}

# ── 2. Bootstrap: BiocManager + remotes ──────────────────────────────────────
miss <- needs(c("BiocManager", "remotes"))
if (length(miss)) { install.packages(miss); clear_locks() }

# ── 3. CRAN dependencies — install in dependency order ───────────────────────
# magick must come before SpatialExperiment (which needs the libmagick++-dev
# system library — install that first: sudo apt-get install -y libmagick++-dev)
cran_layers <- list(
  c("magick", "glmnet", "igraph"),            # level 1 — no inter-deps
  c("DiagrammeR", "ComICS"),                   # level 2 — need igraph/glmnet
  c("purrr", "dplyr", "magrittr", "readr", "readxl",
    "tibble", "rlang", "stringr", "MASS", "matrixStats",
    "data.tree", "limSolve", "e1071", "testit") # level 3 — immunedeconv core
)
for (layer in cran_layers) {
  miss <- needs(layer)
  if (length(miss)) {
    cat("Installing CRAN:", paste(miss, collapse = ", "), "\n")
    install.packages(miss, dependencies = TRUE, Ncpus = parallel::detectCores())
    clear_locks()
  }
}

# ── 4. Bioconductor dependencies ─────────────────────────────────────────────
# Order matters:
#   magick (CRAN, step 3) → SpatialExperiment → GSVA → ConsensusTME (GitHub, step 5)
bioc_deps <- c(
  "BiocParallel",      # parallel evaluation
  "preprocessCore",    # quantile normalisation (EPIC)
  "Biobase",           # ExpressionSet
  "biomaRt",           # gene ID mapping
  "sva",               # surrogate variable correction
  "quantiseqr",        # quanTIseq R implementation
  "SpatialExperiment", # required by GSVA >= 1.52 (Bioc 3.19+)
  "GSVA"               # required by ConsensusTME
)
miss <- needs(bioc_deps)
if (length(miss)) {
  cat("Installing Bioconductor:", paste(miss, collapse = ", "), "\n")
  BiocManager::install(miss, ask = FALSE, update = FALSE,
                       Ncpus = parallel::detectCores())
  clear_locks()
}

# ── 5. GitHub-only dependencies ───────────────────────────────────────────────
# Always force-install to replace any broken partial installs from prior runs.
for (repo in c(
  "GfellerLab/EPIC",
  "dviraran/xCell",
  "grst/MCPcounter",
  "grst/mMCPcounter",
  "cansysbio/ConsensusTME"   # needs GSVA pre-installed (step 4)
)) {
  gh_install(repo)
}

# ── 6. immunedeconv (GitHub only — BiocManager::install will NOT work) ────────
cat("Installing immunedeconv from GitHub...\n")
remotes::install_github("omnideconv/immunedeconv", upgrade = "never", force = TRUE)
clear_locks()

# ── 7. Verify ─────────────────────────────────────────────────────────────────
if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  stop("immunedeconv failed to install. Review errors above.")
}
library(immunedeconv)
cat("\nAvailable deconvolution methods:\n")
print(deconvolution_methods)
cat("\nAvailable TIMER cancer types:\n")
print(timer_available_cancers)
cat("\nSetup complete.\n")
