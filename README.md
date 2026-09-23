# OSA → Heart Failure: Adiposity-Shared Genetic Architecture

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22916622.svg)](https://doi.org/10.5281/zenodo.22916622)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Code repository:** <https://github.com/yl543581-create/osa-heart-failure-mr>
**Archived release (cite this):** <https://doi.org/10.5281/zenodo.22916622>
**Analysis plan:** deposited at OSF (see `docs/OSF_ANALYSIS_PLAN.md`)

Analysis code for a two-sample and multivariable Mendelian randomisation (MR) study
examining whether the genetic association between obstructive sleep apnoea (OSA)
and heart failure is independent of adiposity.

## Key finding

Genetically predicted OSA is robustly associated with heart failure, but this
association **primarily reflects shared genetic architecture with adiposity**
rather than a distinct causal pathway:

| Analysis | Result |
|---|---|
| Univariable MR (MVP → HERMES) | OR 1.279 (1.187–1.378), p = 1.0e-10 |
| Independent replication (meta → FinnGen) | OR 1.345 (1.103–1.641), p = 0.0035 |
| **MVMR adjusted for BMI** | **OR 1.158 (0.904–1.483), p = 0.246** |
| Conditional F for OSA given BMI | **40.9** (adequate power) |
| PheWAS: effect on obesity | **β = +1.95, p = 3.0e-8** (largest of 17 traits) |
| Colocalisation: strongest shared locus | **FTO, PP.H4 = 0.999** |
| Steiger directionality | 10/10 instruments correct direction |

A BMI-stratified analysis partially qualifies the mediation interpretation:
the association persists within the overweight stratum (OR 1.277, p = 0.005),
suggesting residual effects at the upper end of the BMI distribution cannot
be excluded. See `docs/` for full interpretation.

## Data sources (all public)

| Role | Dataset | Accession | N |
|---|---|---|---|
| Exposure (primary) | MVP OSA, European | `GCST90475824` | 152,031 cases / 278,027 controls |
| Exposure (meta + replication) | FinnGen R13 OSA | `G6_SLEEPAPNO_INCLAVO` | 69,677 cases / 430,509 controls |
| Outcome (primary) | HERMES heart failure | `GCST009541` | 47,309 / 930,014 |
| Outcome (replication) | FinnGen R13 heart failure | `I9_HEARTFAIL` | 41,591 / 458,595 |
| Outcome (stratified) | FinnGen HF + BMI≥25 | `I9_HEARTFAIL_AND_OVERWEIGHT` | 25,129 / 204,333 |
| Outcome (stratified) | FinnGen HF + CHD | `I9_HEARTFAIL_AND_CHD` | 26,459 / 409,472 |
| Confounder (MVMR) | BMI | `ieu-b-40` | 501 instruments |

Download URLs:
- GWAS Catalog: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/`
- FinnGen R13: `https://storage.googleapis.com/finngen-public-data-r13/summary_stats/`
- OpenGWAS API: `https://api.opengwas.io/`

## Repository structure

```
.
├── config.R              # portable path resolution
├── _bootstrap.R          # sourced by every script
├── README.md
├── LICENSE
├── CITATION.cff
├── docs/                 # analysis reports and interpretation
├── scripts/              # numbered pipeline (00 → 46)
├── data/                 # NOT in git: downloaded GWAS summary statistics
├── results/              # small result tables
└── logs/                 # run logs
```

## Requirements

- **R ≥ 4.1** (uses the native pipe `|>`); developed on R 4.4.1
- R packages: `TwoSampleMR` (0.7.9), `ieugwasr` (1.1.0), `coloc` (5.2.3),
  `data.table` (1.18.2.1), `MRPRESSO` (1.0), `MendelianRandomization` (0.10.0),
  `RadialMR` (1.2.4)

```r
install.packages(c("data.table","coloc","MRPRESSO","MendelianRandomization",
                   "RadialMR","remotes"))
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/ieugwasr")
```

**OpenGWAS credentials.** Several steps query the OpenGWAS API and require a
token. Set it before running:

```r
Sys.setenv(OPENGWAS_JWT = "<your-token>")
```

Get a token at <https://api.opengwas.io/>.

### Packages that could NOT be installed in the development environment

- **`MVMR`** — unavailable (r-universe returned HTTP 403). MVMR-IVW and the
  Sanderson (2019) conditional F are therefore **implemented directly** in
  `scripts/11_mvmr_final.R` and `scripts/22_mvmr_meta.R`.
- **`MRlap`** — requires `GenomicSEM` → `lavaan` → `sfsmisc`, which need a
  source build toolchain; `make` was absent. Sample-overlap correction is
  therefore not implemented. See `docs/` for the mitigating argument.

## Running the pipeline

All scripts resolve paths relative to their own location, so they can be run
from any working directory.

```bash
Rscript scripts/01_search_opengwas.R        # dataset discovery, instrument counts
Rscript scripts/02_probe_and_main_mr.R      # first-pass MR (FinnGen instrument)
Rscript scripts/04_download_and_extract_mvp_osa.R   # ~548 MB download
Rscript scripts/08_extract_mvp_v2.R         # extract genome-wide significant SNPs
Rscript scripts/09_mvp_instruments_and_mr.R # MVP instruments + main MR
Rscript scripts/13_clump_proper.R           # LD clumping (r2 < 0.001)
Rscript scripts/14_definitive_analysis.R    # definitive main + MVMR analysis
Rscript scripts/17_download_finngen_v2.R    # FinnGen OSA + HF (~1.5 GB)
Rscript scripts/21_meta_final.R             # MVP + FinnGen meta-analysis
Rscript scripts/22_mvmr_meta.R              # MVMR with meta-analysed instrument
Rscript scripts/24_finngen_replication_fast.R  # independent replication
Rscript scripts/31_coloc_local.R            # colocalisation (local data)
Rscript scripts/32_steiger_local.R          # Steiger directionality
Rscript scripts/35_phewas_v2.R              # instrument PheWAS
Rscript scripts/42_overweight_stratified_mr.R  # BMI-stratified endpoints
Rscript scripts/43_collect_metadata.R       # software versions, counts
```

**Scripts 34 and 35 are superseded pairs** — `35` is the corrected version of
`34`; run `35`. Likewise `19/20/21` (run `21`), `25/27/31` (run `31`), and
`26/28/32` (run `32`). Earlier numbers are retained to document the
troubleshooting path.

Total download volume is approximately **5.4 GB**.

## Documents

| File | Content |
|---|---|
| `docs/PAPER_SKELETON.md` | Manuscript outline, figures, discussion points |
| `docs/STROBE_MR_checklist.md` | STROBE-MR reporting checklist |
| `docs/FINAL_CONCLUSION.md` | Consolidated results and evidence chain |
| `docs/OSF_ANALYSIS_PLAN.md` | Analysis plan, as deposited at OSF |

## Known limitations

1. European ancestry only
2. Instruments are enriched at adiposity loci (this is itself the finding)
3. The MVP OSA file's `standard_error` column is entirely `#NA`; standard errors
   were derived from p-values (`se = |beta| / qnorm(p/2)`), because the reported
   confidence intervals are not self-consistent with the p-values
4. The meta-analysed instrument comprises only 10 independent SNPs
5. `MRlap` sample-overlap correction could not be run (see above)
6. Colocalisation was performed within a single cohort (FinnGen vs FinnGen)
7. The study was not prospectively registered; the analysis plan was deposited
   retrospectively at OSF (`docs/OSF_ANALYSIS_PLAN.md`)

## Citation

If you use this code, please cite the archived release:

```bibtex
@software{yang2026osa,
  author    = {Yang, Lujing and Yang, Xiaona},
  title     = {OSA and heart failure: adiposity-shared genetic architecture},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22916622},
  url       = {https://doi.org/10.5281/zenodo.22916622}
}
```

## License

MIT — see `LICENSE`.
