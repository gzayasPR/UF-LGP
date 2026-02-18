# Skim-Seek VCF Merging Pipeline

Merge per-sample Neogen Skim-Seek VCF files with confidence-aware genotype filtering.

## Overview

This pipeline merges individual Skim-Seek VCF files into a multi-sample VCF, with options to handle low-confidence genotype calls. The pipeline handles:
- PASS and LOWCONF filter records
- Genotype probability (GP) based filtering
- Sample renaming and ID standardization
- Parallel processing for large sample sets
- Multi-contig merging optimized for cattle genome

## Pipeline Entry Point

**Script:** `bin/pipelines/skim_seek.sh`

## Key Features

### Confidence-Aware Filtering

**Two modes controlled by `KEEP_LOWCONF`:**

1. **`KEEP_LOWCONF=0` (Recommended):** Set low-confidence genotypes to missing
   - Low-confidence defined as: `FILTER="LOWCONF"` OR `MAX(FMT/GP) < GP_MIN`
   - Maintains data quality by removing uncertain calls
   
2. **`KEEP_LOWCONF=1`:** Keep all genotypes as called
   - Useful for downstream filtering or maximum data retention

## Required Configuration

Edit these variables in `skim_seek.sh`:

```bash
proj_env="/path/to/bin/project_env.sh"
NAME="Skim_Seek_2026"
meta_csv="/path/to/metadata/SSK_IDs_2026.csv"

# Confidence filtering parameters
KEEP_LOWCONF=0        # 0=mask low-confidence, 1=keep all
GP_MIN=0.95           # Minimum genotype probability (used if KEEP_LOWCONF=0)

# SLURM resources
CPUS=30               # Parallel workers
mem_CPU=12Gb
max_time=120:00:00
```

## Input Files

### Sample Metadata CSV (`meta_csv`)

**Required columns:**
- `path` - Full path to per-sample VCF file
- `SkimSeek_ID` - Original sample ID in VCF
- `ProjectID` - Project identifier (optional)
- `IID` - Individual ID (optional)
- `NewID` - Final sample name in merged VCF

**Format:**
```csv
path,SkimSeek_ID,ProjectID,IID,NewID
/blue/data/sample1.vcf.gz,SSK001,Project1,Animal001,UF_Animal001
/blue/data/sample2.vcf.gz,SSK002,Project1,Animal002,UF_Animal002
```

### Per-Sample VCF Requirements

- **Format:** VCF 4.2+, bgzipped (.vcf.gz)
- **Index:** Tabix index (.tbi) preferred (auto-created if missing)
- **Samples:** Exactly 1 sample per file
- **Contigs:** Standard cattle chromosomes (1-29, X, MT, plus NKLS contigs)

**Expected FORMAT fields:**
```
GT  - Genotype (required)
GP  - Genotype probabilities (used for filtering)
AC  - Allele count (stripped during processing)
RC  - Read count (stripped during processing)
DS  - Dosage (stripped during processing)
```

## Processing Workflow

### Stage 1: Metadata Parsing

1. Read CSV and validate columns
2. Check for duplicate `SkimSeek_ID` or `NewID`
3. Verify all VCF files exist
4. Build processing lists

### Stage 2: Per-Sample Cleaning (Parallel)

For each sample:

1. **Filter records:**
   ```bash
   bcftools view -U  # Drop missing GT
   bcftools view -i 'FILTER="PASS" || FILTER="LOWCONF"'  # Keep both
   ```

2. **Optional: Mask low-confidence genotypes**
   
   If `KEEP_LOWCONF=0`:
   ```bash
   bcftools +setGT -- -t q -n . -i 'FILTER="LOWCONF" || MAX(FMT/GP)<GP_MIN'
   ```
   This sets GT to `./. ` for low-confidence calls.

3. **Strip to GT-only:**
   ```bash
   bcftools annotate -x INFO,FORMAT/RC,FORMAT/AC,FORMAT/GP,FORMAT/DS
   ```
   Removes all INFO and extra FORMAT fields (keeps only GT).

4. **Rename sample:**
   ```bash
   bcftools reheader -s [SkimSeek_ID -> SkimSeek_ID]
   ```
   Forces sample name to exact `SkimSeek_ID` from metadata.

**Output:** `clean_vcfs/[SkimSeek_ID].clean.vcf.gz`

### Stage 3: Multi-Contig Merging (Parallel)

1. **Build contig order:**
   - Autosomes: 1-29 (sorted numerically)
   - Sex chromosomes: X
   - Mitochondrial: MT
   - Other: NKLS and non-standard contigs

2. **Merge per contig:**
   ```bash
   bcftools merge --force-samples -r [CONTIG] --file-list [cleaned_vcfs]
   ```

3. **Concatenate contigs:**
   ```bash
   bcftools concat -a --threads [N] -f [parts_list]
   ```

### Stage 4: Sample Renaming

Map `SkimSeek_ID` → `NewID` using metadata:

```bash
bcftools reheader -s [samples_map]
```

### Stage 5: Variant ID Assignment

Set variant IDs to `CHR:POS:REF:ALT`:

```bash
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT'
```

## Output Files

### Final Merged VCF

**Location:** `results/[NAME]/[NAME]_merged_output.vcf.gz`

- Multi-sample VCF with all samples
- GT-only FORMAT (all other fields removed)
- Indexed with tabix
- Variant IDs: `CHR:POS:REF:ALT`
- Sample names: `NewID` from metadata

### Metadata Files

```
results/[NAME]/
├── samples_2026.tsv              # path, SkimSeek_ID, NewID mapping
├── vcf_paths.list                # Original VCF paths
├── vcf_to_merge.cleaned.list     # Cleaned VCF paths
├── _current_samples.txt          # Sample order before renaming
├── _new_samples.txt              # Sample order after renaming
└── [NAME]_stats.txt              # bcftools stats output
```

### Logs

**Location:** `bash_out/Skim_Seek/`

- Main job log: `Skim_Seek_2026_[JOBID].log`
- Per-sample cleaning logs: `logs/[SkimSeek_ID].clean.log`

## Quality Control

### GP Score Distribution

If using `KEEP_LOWCONF=0`, check how many genotypes are masked:

```bash
# Before filtering
bcftools stats [original.vcf.gz] | grep "number of samples"

# After filtering  
bcftools stats [merged_output.vcf.gz] | grep "number of samples"
bcftools stats [merged_output.vcf.gz] | grep "number of SNPs"
```

### Missingness Rates

```bash
# Per-sample missingness
bcftools query -f '[%SAMPLE\t%GT\n]' [merged_output.vcf.gz] | \
  awk '$2=="./." {miss[$1]++} {total[$1]++} END {for (s in total) print s, miss[s]/total[s]}'

# Per-variant missingness
bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' [merged_output.vcf.gz] | \
  awk '{miss=0; for(i=4;i<=NF;i++) if($i=="./.") miss++; print $1,$2,$3,miss/(NF-3)}'
```

## Parameter Tuning

### GP_MIN Threshold

**Conservative (GP_MIN=0.95):**
- Keeps only high-confidence genotypes
- Higher missingness but fewer errors
- Recommended for genomic prediction

**Moderate (GP_MIN=0.90):**
- Balance between data retention and quality
- Standard for most analyses

**Permissive (GP_MIN=0.80):**
- Retains more data
- Use with caution for critical analyses

### KEEP_LOWCONF Decision

| Analysis Type | Recommended Setting | Rationale |
|--------------|---------------------|-----------|
| Genomic prediction | KEEP_LOWCONF=0, GP_MIN=0.95 | Minimize genotype errors |
| Population structure | KEEP_LOWCONF=0, GP_MIN=0.90 | Balance quality and coverage |
| Rare variant detection | KEEP_LOWCONF=1 | Maximize data, filter downstream |
| Imputation reference | KEEP_LOWCONF=0, GP_MIN=0.95 | High-quality reference critical |

## SLURM Resource Scaling

| Sample Count | CPUs | Memory/CPU | Runtime (est) |
|-------------|------|-----------|---------------|
| <100 | 12 | 8GB | 2-4 hours |
| 100-500 | 20 | 12GB | 4-8 hours |
| 500-1000 | 30 | 12GB | 8-16 hours |
| >1000 | 40 | 16GB | 16-48 hours |

## Common Issues

### Issue 1: Duplicate Sample Names

**Error:** `Duplicate SkimSeek_ID in metadata`

**Solution:** Check metadata CSV for duplicate entries in `SkimSeek_ID` or `NewID` columns

### Issue 2: Missing bcftools setGT Plugin

**Error:** `bcftools plugin 'setGT' not available`

**Solution:**
- Load bcftools with plugins: `ml bcftools/1.22` (or compiled with plugins)
- OR set `KEEP_LOWCONF=1` to skip filtering

### Issue 3: Contig Name Mismatches

**Error:** Contigs missing after merge

**Solution:** Check contig names in VCF headers - ensure consistent naming (e.g., "1" vs "chr1")

### Issue 4: Memory Exhaustion

**Error:** `Killed` or memory errors

**Solution:**
- Increase `--mem-per-cpu` in SLURM header
- Reduce `CPUS` to lower total memory
- Process in chromosome batches

## Usage Example

```bash
# 1. Create metadata CSV
cat > metadata/SSK_IDs_2026.csv << EOF
path,SkimSeek_ID,ProjectID,IID,NewID
/blue/data/SSK001.vcf.gz,SSK001,P1,A001,UF_A001
/blue/data/SSK002.vcf.gz,SSK002,P1,A002,UF_A002
EOF

# 2. Edit configuration in skim_seek.sh
vim bin/pipelines/skim_seek.sh
# Set KEEP_LOWCONF=0, GP_MIN=0.95

# 3. Submit job
sbatch bin/pipelines/skim_seek.sh

# 4. Monitor progress
tail -f bash_out/Skim_Seek/Skim_Seek_2026_*.log

# 5. Check output
bcftools stats results/Skim_Seek_2026/Skim_Seek_2026_merged_output.vcf.gz
```

## Dependencies

- **bcftools 1.22+** - VCF manipulation (must include +setGT plugin if KEEP_LOWCONF=0)
- **htslib** - tabix, bgzip
- **GNU parallel** (optional) - Faster parallel processing
- **samtools** (optional) - Additional VCF utilities

## Performance Optimization

**Parallel Processing:**
- Use GNU parallel if available (faster than xargs)
- Scale `CPUS` based on sample count
- Per-contig merging parallelizes naturally

**Disk I/O:**
- Use local scratch space if available: `scratch_root="${TMPDIR}"`
- Keep intermediate files on fast storage

**Memory:**
- Scales with samples × SNPs in memory during merge
- Typical: 1000 samples × 500K SNPs ≈ 10-15GB

## Next Steps

After successful merging:
- Convert to PLINK for downstream analysis
- Merge with Illumina 250K data (see `cross_platform.md`)
- Perform QC filtering and imputation
- Run population structure or association analyses
