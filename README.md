# UF-LGP: University of Florida - Livestock Genomics Pipelines

A comprehensive collection of genomics analysis pipelines for livestock research, developed during PhD and postdoctoral work at the University of Florida. These pipelines provide end-to-end workflows for processing and analyzing cattle genotype data from multiple platforms.

## Overview

UF-LGP provides production-ready pipelines for:
- **Multi-platform genotype processing**: Illumina SNP arrays and Neogen Skim-Seek sequencing
- **Breed composition analysis**: ADMIXTURE-based ancestry inference with visualization
- **Parentage verification**: KING-based relatedness and parent-offspring identification
- **Population genomics**: PCA, relatedness filtering, and quality control workflows

All pipelines are designed for SLURM-managed HPC environments and handle large-scale cattle genomics datasets.

---

## Repository Structure

```
UF-LGP/
├── bin/
│   ├── pipelines/          # Main pipeline entry points
│   ├── scripts/            # Core processing scripts
│   └── project_env.sh      # Environment configuration
├── metadata/
│   ├── IDs/                # Sample ID lists for analysis
│   └── Illumina_to_Plink/  # Reference files for array processing
├── docs/                   # Documentation (in development)
├── results/                # Pipeline outputs (not tracked)
└── bash_out/               # SLURM logs (not tracked)
```

---

## Pipelines

### 1. Illumina SNP Array Processing (`illumina_2_plink.sh`)

Converts Illumina Bovine 250K FinalReport files to PLINK format with proper genome coordinates.

**Features:**
- Batch processing of multiple FinalReport files
- Automatic format validation and error detection
- ID standardization and coordinate mapping to ARS-UCD1.2
- Handles multi-date genotyping batches

**Usage:**
```bash
sbatch bin/pipelines/illumina_2_plink.sh
```

**Inputs:**
- `metadata/250K_finalreports.csv`: Paths to Illumina FinalReport files
- `metadata/Illumina_to_Plink/`: SNP maps, chromosome updates, allele reference files

**Outputs:**
- `results/250K_YYYY/UFID_250K_YYYY_Illumina.ID_ARS.{bed,bim,fam}`: PLINK binary files

---

### 2. Neogen Skim-Seek VCF Merging (`skim_seek.sh`)

Merges per-sample Skim-Seek VCF files with confidence-aware genotype filtering.

**Features:**
- Handles PASS and LOWCONF filter records
- Optional low-confidence genotype masking (based on GP scores)
- Parallel processing with automatic sample renaming
- Multi-contig merging optimized for cattle genome

**Usage:**
```bash
sbatch bin/pipelines/skim_seek.sh
```

**Configuration:**
- `KEEP_LOWCONF=0`: Mask low-confidence genotypes as missing
- `GP_MIN=0.95`: Minimum genotype probability threshold

**Inputs:**
- `metadata/SSK_IDs_2026.csv`: Sample metadata with VCF paths and IDs

**Outputs:**
- `results/Skim_Seek_YYYY/[NAME]_merged_output.vcf.gz`: Multi-sample VCF
- Sample statistics and QC reports

---

### 3. Cross-Platform Genotype Integration (`SSK_250K.sh`)

Merges Skim-Seek and Illumina 250K datasets on common SNP sites.

**Features:**
- Identifies overlapping SNPs between platforms
- Normalizes and harmonizes allele coding
- Maintains sample union with site intersection
- Quality-filtered common variant set

**Usage:**
```bash
sbatch bin/pipelines/SSK_250K.sh
```

**Outputs:**
- `results/[NAME]/merge_SSK_250K/SkimSeek_250K.commonSites.merged.vcf.gz`

---

### 4. Breed Composition Analysis (`breed_composition.sh`)

Two-stage supervised ADMIXTURE analysis for ancestry inference.

**Stage 1: Population Preparation (`gbc_setup_unrelated.sh`)**
- QC filtering (MAF, missingness, HWE)
- LD pruning for independent SNP set
- PC-AiR based unrelated individual selection
- SNP set harmonization between reference and target populations

**Stage 2: ADMIXTURE Analysis (`gbc_admixture.sh`)**
- Supervised K=2 ancestry estimation
- Reference allele frequency estimation from purebred populations
- Generates Q-plots and PCA visualizations with ancestry overlay

**Usage:**
```bash
sbatch bin/pipelines/breed_composition.sh
```

**Inputs:**
- Reference population IDs (purebred animals)
- Mixed population IDs (animals to analyze)

**Outputs:**
- `results/[NAME]/ADMIXTURE/Supervised_Unrelated/[NAME].combined.ancestry.csv`
- `results/[NAME]/Figures/admixture.png`: Stacked bar plot
- `results/[NAME]/Figures/PCA_Admixture.png`: PCA with ancestry proportions

---

### 5. Parentage Verification (`check_parentage.sh`)

KING-based relatedness inference for parent-offspring identification.

**Features:**
- Quality-controlled SNP pruning for relatedness estimation
- Batch processing of query individuals
- Automatic parent-offspring pair extraction (InfType='PO')
- Summary reports with kinship coefficients and IBS0 statistics

**Usage:**
```bash
sbatch bin/pipelines/check_parentage.sh
```

**Inputs:**
- `metadata/check.parentage.ID`: List of individuals to verify
- Reference population for comparison

**Outputs:**
- `results/Parentage_YYYY/parent_offspring_summary.csv`: PO pairs with statistics
- Per-individual `.kin` files with all relationships
- KING kinship matrix for downstream analysis

**Interpretation:**
- Parent-offspring pairs: Kinship ≈ 0.25, IBS0 ≈ 0
- Full-siblings: Kinship ≈ 0.25, IBS0 > 0
- Second-degree: Kinship ≈ 0.125

---

## Core Utilities

### Shell Functions (`scripts/shell.functions.sh`)

**`detect_genotype_format()`**
- Auto-detects input format: VCF (.vcf, .vcf.gz, .bcf), PLINK1 (.bed/.ped), PLINK2 (.pgen)
- Returns appropriate PLINK flag for seamless format handling
- Used across all pipelines for flexible input

### R Scripts

**`Selection_of_Unrelated_script.R`**
- PC-AiR (Principal Components Analysis in Related samples)
- KING kinship-based relatedness filtering
- Generates unrelated sample lists for population structure analysis

**`Q_PLOT.R`**
- ADMIXTURE ancestry bar plots
- Highlights unrelated individuals used in supervised analysis
- Customizable color schemes for breed groups

**`PCA_plots.R`**
- PCA visualization with ancestry proportions
- Color-coded by breed composition
- Interactive legends and sample labels

**`combine_Q_FAM.R`**
- Merges ADMIXTURE Q-matrix with PLINK FAM files
- Creates analysis-ready ancestry tables

---

## Setup and Configuration

### Initial Setup

```bash
# Clone repository
git clone https://github.com/gzayasPR/UF-LGP.git
cd UF-LGP

# Run setup script to create directory structure
bash Setup.sh

# This creates:
#   - bin/project_env.sh (environment variables)
#   - metadata/, results/, docs/, bash_out/ directories
```

### Environment Configuration

The `bin/project_env.sh` file is auto-generated and contains:
```bash
export proj_dir="/path/to/UF-LGP"
export my_bin="${proj_dir}/bin"
export my_metadata="${proj_dir}/metadata"
export my_results="${proj_dir}/results"
export my_docs="${proj_dir}/docs"
export my_bash="${proj_dir}/bash_out"
```

All pipelines source this file for consistent path management.

---

## Requirements

### Software Dependencies

- **PLINK** 1.90b3+ (genotype QC and format conversion)
- **PLINK2** (modern format support)
- **ADMIXTURE** 1.3.0+ (ancestry inference)
- **KING** 2.3.0+ (relatedness estimation)
- **bcftools** 1.22+ (VCF manipulation)
- **htslib** (tabix, bgzip)
- **R** 3.6+ with packages:
  - `SNPRelate` (GDS format, PCA)
  - `GENESIS` (PC-AiR)
  - `gdsfmt`
  - `ggplot2`, `dplyr`, `tidyr`

### System Requirements

- SLURM workload manager
- 8-32 GB RAM per job (depends on dataset size)
- Multi-core CPU for parallel processing

---

## Input Data Format

### Sample ID Lists
Tab or space-delimited, 2 columns (FID, IID):
```
Family1  Animal001
Family1  Animal002
```

### Illumina Metadata CSV
Header: `full_path`
```
full_path
/path/to/Univ_of_Florida_BOVF250V1_20240101_FinalReport.txt
```

### Skim-Seek Metadata CSV
Header: `path,SkimSeek_ID,ProjectID,IID,NewID`
```
path,SkimSeek_ID,ProjectID,IID,NewID
/path/to/sample1.vcf.gz,SSK001,Project1,Animal001,UF_Animal001
```

---

## Quality Control Parameters

Default QC thresholds (modifiable via pipeline arguments):

- **Individual missingness** (`--mind`): 0.1 (10%)
- **SNP missingness** (`--geno`): 0.01 (1%)
- **Minor allele frequency** (`--maf`): 0.01-0.05
- **Hardy-Weinberg equilibrium**: p < 0.00001
- **LD pruning**: 50 SNP windows, r² < 0.2-0.5

---

## Output Files

### Genotype Files
- **PLINK binary**: `.bed/.bim/.fam` (SNP array data)
- **VCF**: `.vcf.gz` + `.tbi` (sequence data)

### Analysis Results
- **Ancestry**: `.Q` (ancestry proportions), `.P` (allele frequencies)
- **Relatedness**: `.kin` (kinship matrix), `.kin0` (pairwise relationships)
- **PCA**: CSV files with PC scores

### Visualizations
- **PNG plots**: Admixture bar plots, PCA scatterplots with ancestry overlay

---

## Example Workflows

### Complete Breed Composition Analysis

```bash
# 1. Process Illumina genotypes
sbatch bin/pipelines/illumina_2_plink.sh

# 2. Run breed composition pipeline
# Edit breed_composition.sh to set:
#   - IN_GENO: path to PLINK fileset
#   - reflist: purebred reference IDs
#   - mixedlist: animals to analyze

sbatch bin/pipelines/breed_composition.sh

# Outputs:
#   results/[NAME]/Figures/admixture.png
#   results/[NAME]/ADMIXTURE/Supervised_Unrelated/[NAME].combined.ancestry.csv
```

### Parentage Verification Workflow

```bash
# 1. Create list of animals to verify (one ID per line)
echo "2190762" > metadata/check.parentage.ID
echo "2210676" >> metadata/check.parentage.ID

# 2. Run parentage pipeline
sbatch bin/pipelines/check_parentage.sh

# 3. Check results
cat results/Parentage_YYYY/parent_offspring_summary.csv
#   Individual,Parent,InfType,Kinship,IBS0
#   2190762,1040084,PO,0.2456,0.0012
```

---

## Troubleshooting

### Common Issues

**1. PLINK format errors**
- Ensure `.bed`, `.bim`, `.fam` files have matching sample/SNP counts
- Use `detect_genotype_format()` function for automatic detection

**2. VCF indexing failures**
- Re-compress with `bcftools view -Oz` (BGZF format required)
- Re-index with `tabix -p vcf`

**3. ADMIXTURE convergence issues**
- Increase LD pruning stringency (lower r² threshold)
- Remove SNPs with high missingness
- Verify reference population labels

**4. Memory errors**
- Increase `--mem-per-cpu` in SLURM header
- Split large VCFs by chromosome before merging

---

## Citation

If you use these pipelines in your research, please cite:

```
Zayas, G. (2026). UF-LGP: University of Florida Livestock Genomics Pipelines.
GitHub repository: https://github.com/gzayasPR/UF-LGP
```

---

## Contact

**Author:** G. Zayas  
**Institution:** University of Florida, Department of Animal Sciences  
**Email:** gzayas97@ufl.edu

For bug reports or questions about pipeline usage, please open an issue on GitHub.

---

## License

This software is provided for academic and research use. Please contact the University of Florida for licensing information regarding commercial use.

---

## Acknowledgments

Developed in the Mateescu Lab at the University of Florida. These pipelines were created to support livestock genomics research and are provided to the university for continued use and development.

**Genome Assembly:** ARS-UCD1.2 (Bos taurus reference)  
**SNP Arrays:** Illumina Bovine 250K, Neogen Skim-Seek  
**Analysis Software:** PLINK, ADMIXTURE, KING, bcftools, R/Bioconductor
