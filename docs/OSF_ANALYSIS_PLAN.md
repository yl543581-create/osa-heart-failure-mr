# Analysis Plan — OSA and Heart Failure Mendelian Randomisation Study

> **Status: RETROSPECTIVE (post hoc) deposit.**
>
> This document records the analysis plan for a study whose data analysis has
> already been completed. It is **not** a prospective pre-registration and must
> not be described as one. It is deposited to document that the analytical
> choices below were specified in writing before the outcome analyses were
> examined, and to allow readers to assess which results were confirmatory
> versus exploratory.
>
> Deposit date: 2026-09-21
> OSF link: (project URL assigned on submission; see the registration page)

---

## 1. Study question

Is the genetic association between obstructive sleep apnoea (OSA) and heart
failure independent of adiposity, and does adiposity explain that association?

## 2. Motivation

Observational studies consistently associate OSA with heart failure, but
Mendelian randomisation (MR) studies disagree:

- **PMID 36480010** (*Sleep* 2023): multivariable MR adjusting for BMI found the
  OSA–heart-failure association persisted (IVW OR 1.13, 1.01–1.27), using
  **5 genome-wide significant SNPs**.
- **PMID 36611115** (*Eur J Prev Cardiol* 2023): the same adjustment strategy
  found **no** association.
- **PMID 40472801** (*EBioMedicine* 2025): a polygenic **score** study (not MR)
  found the heart-failure signal disappeared once BMI genetic contributions were
  removed from the OSA score.
- **PMID 39288744** (*Cardiology* 2025): two-step MR identified obesity, glucose
  and depression as mediators, without BMI-adjusted MVMR.

The contradiction has not been adjudicated. A plausible explanation is that the
earlier instruments were too weak to resolve the question.

## 3. Hypotheses

- **H1**: The genetic association between OSA and heart failure persists after
  adjustment for BMI (MVMR).
- **H2**: The OSA instrument explains more variance in OSA than in heart failure
  (Steiger directionality).
- **H3**: The OSA instrument is functionally an adiposity instrument
  (PheWAS and colocalisation).
- **H4**: Stratifying the heart-failure outcome by adiposity status changes the
  effect size.

## 4. Data sources (fixed before analysis)

| Role | Dataset | Accession |
|---|---|---|
| Exposure (primary) | MVP OSA, European ancestry | `GCST90475824` |
| Exposure (meta + replication) | FinnGen R13 sleep apnoea incl. outpatient | `G6_SLEEPAPNO_INCLAVO` |
| Outcome (primary) | HERMES heart failure | `GCST009541` |
| Outcome (replication) | FinnGen R13 heart failure | `I9_HEARTFAIL` |
| Outcome (stratified) | FinnGen R13 heart failure + BMI≥25 | `I9_HEARTFAIL_AND_OVERWEIGHT` |
| Outcome (stratified comparator) | FinnGen R13 heart failure + CHD | `I9_HEARTFAIL_AND_CHD` |
| Confounder for MVMR | Body mass index | `ieu-b-40` |

**Population**: European ancestry only.

## 5. Pre-specified instrument selection

| Parameter | Value |
|---|---|
| Significance threshold | **P < 5 × 10⁻⁸** (not relaxed) |
| LD clumping | **r² < 0.001** |
| Clumping window | **10,000 kb** |
| Reference panel | 1000 Genomes EUR (via OpenGWAS API) |
| Minimum instrument strength | **F > 10** |

**Sensitivity variants of the instrument set** (all pre-specified):
1. MVP only
2. MVP + FinnGen meta-analysed (inverse-variance fixed effects)
3. Meta-analysed set clumped to independence

## 6. Pre-specified statistical analyses

### Primary
- Two-sample MR: **inverse-variance weighted (IVW)**
- Report also: MR-Egger, weighted median, weighted mode
- Effect scale: **odds ratio** per 1 log-odds of genetic liability

### Secondary
- **Multivariable MR** adjusting for BMI, with **conditional F statistics**
  (Sanderson et al. 2019)
- **Replication** in FinnGen heart failure
- **Colocalisation** (`coloc.abf`) at each instrument locus, ±500 kb window
- **Steiger directionality** on the liability scale
- **Instrument PheWAS** across pre-specified traits
- **BMI-stratified outcome** comparison

### Pre-specified sensitivity analyses
Cochran Q, MR-Egger intercept, MR-PRESSO, leave-one-out, FTO-region exclusion.

## 7. Pre-specified thresholds

| Quantity | Threshold |
|---|---|
| Primary significance | α = 0.05 |
| Multiple-testing (5 primary tests) | Bonferroni α = 0.01 |
| PheWAS | Benjamini–Hochberg FDR < 0.05 |
| Colocalisation | **PP.H4 ≥ 0.8** considered strong evidence of a shared causal variant; **PP.H3 ≥ 0.8** considered strong evidence of distinct causal variants |
| Instrument strength | F > 10 |
| Conditional F (MVMR) | > 10 considered adequate |

## 8. Pre-specified interpretation rules

- If the MVMR direct effect is null **and** conditional F > 10 → conclude the
  association is not independent of adiposity, with adequate power.
- If the MVMR direct effect is null **but** conditional F < 10 → conclude the
  analysis is underpowered and no conclusion can be drawn.
- Null results in PheWAS for traits unrelated to adiposity are treated as
  evidence of instrument specificity.

## 9. Deviations and post hoc analyses (declared)

The following were **not** pre-specified and should be treated as exploratory:

1. **Colocalisation** was performed within a single cohort (FinnGen OSA vs
   FinnGen heart failure) rather than across cohorts, because the OpenGWAS
   regional query interface and missing effect-allele frequencies in the
   heart-failure outcome prevented a cross-cohort implementation.
2. **MVMR implementation**: the `MVMR` R package could not be installed
   (r-universe HTTP 403), so MVMR-IVW was implemented directly as weighted least
   squares with Sanderson (2019) conditional F statistics.
3. **Meta-analysis** of MVP and FinnGen OSA was added during the analysis when
   the MVP-only instrument gave a conditional F of 1.3 (inadequate). The
   meta-analysed instrument raised this to 40.9.
4. **Standard errors** in `GCST90475824` are missing (all `#NA`); SE were derived
   from p-values as `se = |beta| / qnorm(p/2)` after establishing that the
   reported confidence intervals are not self-consistent with the p-values.
5. **MRlap** sample-overlap correction could not be performed (build toolchain
   unavailable).

## 10. Software

R 4.4.1; TwoSampleMR 0.7.9; ieugwasr 1.1.0; coloc 5.2.3; data.table 1.18.2.1;
MRPRESSO 1.0; MendelianRandomization 0.10.0; RadialMR 1.2.4.
Full session information: `results/sessionInfo.txt`.

## 11. Reporting

Results are reported following **STROBE-MR** (Skrivankova et al., *JAMA* 2021).
Checklist: `docs/STROBE_MR_checklist.md`.

## 12. Registration statement for the manuscript

Suggested wording for the limitations section:

> "This study was not prospectively registered. The analysis plan, including the
> instrument selection thresholds, statistical methods, significance thresholds
> and interpretation rules, was specified in writing before the outcome analyses
> were conducted and is available at [OSF link]. Analyses that were not part of
> that plan are identified as exploratory in the Methods."
