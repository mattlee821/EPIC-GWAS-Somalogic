# EPIC Somalogic GWAS Pipeline

**Version:** March 2026  

This pipeline runs protein GWAS across multiple case‑control studies, stratified by analysis group (combined, cases, controls). It performs QC, GWAS, post‑GWAS QC/plots, and optional meta‑analysis (within‑study cases+controls and cross‑study meta).

---

## Quick Start

```bash
# Example with comma-separated list
bash src/004_run.sh \
  --profile slurm \
  --out "$ANALYSIS_OUT" \
  --studies "Neuro_01,Brea_01_Erneg" \
  --proteins "10620-21,18878-15,4721-54,8323-163" \
  --pheno "data/phenotype/phenofile-test.txt" \
  --covariates "Sex,Age_Recr" \
  --meta-group true \
  --meta-study "cases,controls,combined,meta"

# Example with file-based list (one per line)
bash src/004_run.sh \
  --proteins "data/phenotype/features.txt" \
  # ... other args
```

---

## Key Inputs

### `pipeline/samplesheet.csv`
Required columns:
- `study_id`
- `plink_bfile` (prefix to .bed/.bim/.fam or .pgen/.pvar/.psam)
- `bgen_dir` (directory with BGENs, if used)
- `sample_file` (sample sheet: `.psam`, `.fam`, or `.sam`)
- `group_column` (e.g. `PHENO` in the sample file)
- `cases_value` (e.g. `2` for case‑control coding)

The pipeline uses `group_column` to split **cases**, **controls**, and **combined** groups.
For the **combined** analysis, the `group_column` (e.g., `PHENO`) is also added as a covariate in the model.

### Global files
Set in `pipeline/params.yaml`, or override at runtime:
- `phenotype_file`
- `covariate_file`

---

## Running the Pipeline (`src/004_run.sh`)

Arguments supported by `src/004_run.sh`:

- `--profile`: Nextflow profile (default: `slurm`)
- `--out`, `-o`, `--outdir`: output directory
- `--studies`: comma‑separated study IDs (subset)
- `--proteins`: subset of traits. Supports:
    - Comma-separated list (e.g. `prot1,prot2`)
    - **Path to a file** (e.g. `protein_list.txt`) with one ID per line
    - **Omit entirely** to run ALL numeric traits in the phenotype file
- `--pheno`: phenotype file path
- `--covar`: covariate file path
- `--covariates`: comma‑separated covariate names (e.g. `Sex,Age_Recr`)
- `--meta-group`: `true|false`, run within‑study meta of cases+controls
- `--meta-study`: comma‑separated list for study‑level meta (any combo of `cases,controls,combined,meta`). Default: `false`

The script always runs Nextflow with `-resume`.

---

## Workflow in `src/`

Primary helper scripts (in order of typical use):

- `000_env.sh`: environment setup (shell helpers)
- `000_tools.sh`: download and install the METAL binary into `tools/`
- `001_data-genetics.sh`: ingest/prepare genotype data paths
- `002_data-phenotype.R`: ingest/prepare phenotype data
- `003_samplesheet.R`: build `pipeline/samplesheet.csv` from metadata
- `004_run.sh`: main pipeline launcher (Nextflow wrapper)
- `005_submit.sh`: SLURM submission wrapper for `004_run.sh`
- `006_clear_meta_cache.sh`: remove Nextflow cache for specific processes (e.g., meta)

Additional:
- `src/test/`: ad‑hoc testing utilities (e.g., meta‑analysis testing)
- `src/archive/`: archived components kept for reference

---

## Pipeline Steps and Filters

Below are the major steps and **all filtering thresholds** used. Thresholds come from `pipeline/params.yaml`.

### 1. Sample QC (PLINK2)
**Script:** `pipeline/bin/sample_qc.sh`  
**Inputs:** genotype data + optional keep list  
**Filters:**
1. **Call rate filter:** `--mind 0.05`  
2. **Heterozygosity outliers:** F‑coefficient within **mean ± 3 SD**  
   - Computed on LD‑pruned markers  
3. **Relatedness filter:** `--king-cutoff 0.0884` (≈ 2nd‑degree relatedness)  

**Output:** `qc_pass_samples.txt`

---

### 2. Variant QC (PLINK2)
**Script:** `pipeline/bin/variant_qc.sh`  
**Filters:**
- `--maf 0.01`
- `--hwe 1e-15`
- `--geno 0.05`

**Output:** `variant_qc_pass.snplist`

---

### 3. LD Pruning
**Script:** `pipeline/bin/ld_pruning.sh`  
**Parameters:**
- `--indep-pairwise 1000kb 1 0.5`  

This creates LD‑pruned variants for REGENIE Step 1 and outputs a pruned PLINK dataset (`step1_input.*`).

---

### 4. PCA (PLINK2)
**Script:** `pipeline/bin/compute_pcs.sh`  
**Settings:** `--pca 10`  

PCs are computed on the **same LD‑pruned data** used for REGENIE Step 1.  
They are merged into the covariate file and used in **both Step 1 and Step 2**.

---

### 5. REGENIE Step 1
**Script:** `pipeline/bin/run_regenie_step1.sh`  
**Key settings:**
- `--step 1`
- `--bsize 1000`
- `--lowmem`
- LOCO predictors generated (`loco_*.loco.gz`)

REGENIE Step 1 expects strictly numeric data. The pipeline automatically:
1. Filters for numeric columns in the phenotype file.
2. **Excludes** columns identified as covariates (`--covariates`).
3. **Excludes** common metadata columns (`STUDY`, `PlateId`, `FID`, `IID`).

**Inputs:** LD‑pruned genotype, covariates, cleaned phenotypes  

---

### 6. REGENIE Step 2
**Script:** `pipeline/bin/run_regenie_step2.sh`  
**Key settings:**
- `--step 2`
- `--bsize 400`
- `--qt`
- `--gz`
- Per‑chromosome and per‑chunk association testing

**Outputs:** raw REGENIE `.regenie.gz` summary stats (per chunk per chromosome)

---

### 7. Post‑GWAS QC + Plots (Python)
**Script:** `pipeline/bin/finalise.py`  
**Functions:**
1. Concatenates all REGENIE outputs for a protein and group
2. Computes **P = 10^(-LOG10P)**
3. Keeps key columns
4. Generates **pre‑QC Manhattan + QQ plots**
5. Filters (QC):
   - INFO score must be at least 0.3.
   - Variants must have sample size at least half of the maximum sample size observed in that group; this removes variants with unusually small effective N.
   - Effect sizes with absolute value greater than 10 are removed to exclude extreme, likely unstable estimates.
   - Standard errors must be positive and not excessively large (kept only if at most 10), which removes invalid or highly unstable regression results.
6. Writes final `gwas.tsv.gz`
7. Generates **post‑QC Manhattan + QQ plots**
8. Writes `metrics.tsv`

**Final GWAS columns (order):**  
`CHR, POS, ID, EA, OA, EAF, BETA, SE, P, N, CHISQ`

---

### 8. Meta‑Analysis (METAL)
**Script:** `pipeline/bin/run_metal.sh`  
**Settings:**
- `SCHEME STDERR`
- `AVERAGEFREQ ON`
- `MINMAXFREQ ON`
- `MARKER ID`
- `ALLELE EA OA`
- `EFFECT BETA`
- `STDERR SE`
- `PVALUE P`
- `WEIGHT N`
- `FREQ EAF`
- **Custom variable to track total N across studies:**
  - `CUSTOMVARIABLE TotalSampleSize`
  - `LABEL TotalSampleSize as N`

**Two meta modes:**
1. **Within‑study meta (`--meta-group true`)**  
   Cases + Controls → `meta/METAL` under each study/protein

2. **Study‑level meta (`--meta-study`)**  
   Any combination of `cases,controls,combined,meta`  
   Output path: `meta/<protein>/<analysis_type>/METAL`

**Meta QC (post‑METAL):**
- `finalise.py` is run on all METAL outputs to generate:
  - `gwas.tsv.gz`
  - `metrics.tsv`
  - `figure.png`
- These are written to a **QC** subdirectory alongside METAL output:
  - Within‑study meta: `<study>/<protein>/meta/QC`
  - Study‑level meta: `meta/<protein>/<analysis_type>/QC`

**Outputs (kept):**
`meta.tbl.gz`, `meta.tbl.info`, `meta.log` (raw METAL)  
`gwas.tsv.gz`, `metrics.tsv`, `figure.png` (QC output)

---

## Output Structure (Summary)

Top‑level:
```
<outdir>/
├── _shared/001_validation
├── <study_id>/
│   ├── _shared/<group>/...            # Steps 2–7 shared outputs
│   └── <protein>/<group>/
│       ├── 008_regenie-step2
│       └── 009_QC
└── meta/
    └── <protein>/<analysis_type>/
        ├── METAL
        └── QC
```

---

## Notes
- `--meta-study` defaults to `false`.  
  You must explicitly set it to run study‑level meta.
- `--meta-group` defaults to `false`.  
  Enable it to create within‑study cases+controls meta.

---

## Scaling to Thousands of Proteins

To run a large-scale analysis:
1. **Provide a list file**: Use `--proteins "/path/to/list.txt"` (one ID per line).
2. **Run All**: Omit `--proteins` entirely. The pipeline will automatically identify and run all 3,000+ numeric traits while ignoring metadata.
3. **Chunking**: By default, Step 2 is chunked (`chunk_size: 100` in `params.yaml`). This splits 3,000 proteins into 30 parallel Slurm jobs for maximum efficiency.
