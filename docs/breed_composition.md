# Breed Composition Analysis Pipeline

Supervised ADMIXTURE analysis for inferring breed ancestry proportions with PCA visualization.

## Overview

This two-stage pipeline estimates ancestry proportions in admixed populations using:
1. **Unrelated individual selection** - PC-AiR to identify unrelated samples
2. **Supervised ADMIXTURE** - K=2 ancestry estimation with reference populations

Outputs include ancestry proportion estimates and visualizations (Q-plots, PCA with ancestry overlay).

## Pipeline Entry Point

**Script:** `bin/pipelines/breed_composition.sh`

This submits two sequential jobs:
- `gbc_setup_unrelated.sh` - Population preparation
- `gbc_admixture.sh` - ADMIXTURE analysis and visualization

## Required Configuration

Edit these variables in `breed_composition.sh`:

```bash
proj_env="/path/to/bin/project_env.sh"
NAME="GBC_February_2026"
IN_GENO="/path/to/genotypes"     # PLINK fileset (bed/bim/fam)
mixedlist="/path/to/admix.ID"    # Animals to analyze
reflist="/path/to/MAB.ID"        # Purebred reference animals

# QC parameters
mind=0.1    # Max individual missingness
maf=0.01    # Minor allele frequency
geno=0.1    # Max SNP missingness
```

## Input Files

### 1. Genotype Data (`IN_GENO`)

**Format:** PLINK binary fileset (`.bed/.bim/.fam`)

**Requirements:**
- Autosomes 1-29 (sex chromosomes excluded)
- Common SNP set across all samples
- Pre-QC recommended (remove low-quality samples/SNPs)

### 2. Reference Population IDs (`reflist`)

**Format:** Two-column space/tab-delimited (FID, IID)

**For UF(Gabe's comment):** I use the MAB + Brahman animals as the reference population. The IDs I currently use are found in the metadata/IDs/MAB.ID. 
```
MAB_BRA	2130309
MAB_BRA	2130326
MAB_BRA	2130336
Brahman	932619
Brahman	932787
Brahman	935186
Brahman	936741
```

**Requirements:**
- Animal's FID and IID must match the ones in the PLINK files.

### 3. Mixed Population IDs (`mixedlist`)

**Format:** Two-column space/tab-delimited (FID, IID)

**For UF(Gabe's comment):** I use the MAB + Brahman but also include our commerical Brangus populations metadata/IDs/admix.ID. 
```
MAB_BRA	    3999288
SeminNutri	1409
SeminNutri	2447
SeminThermo	500012
SeminThermo	500013
ThermoWes	  170006
ThermoWes	  170009
Brahman   	936741
```

**Requirements:**
- Animals with unknown or mixed ancestry
- Can overlap with `reflist` (will be in both analyses)

## Processing Workflow

### Stage 1: Population Preparation (`gbc_setup_unrelated.sh`)


Based on previous analysis, using an unrelated subset of the MAB is the best for PCA and ADMIXTURE/breed composition analysis. 
Too implement this we need to identify the unrelated subset, which is this current step. For here we do quality control on both the ref and admix population list.

#### Step 1: Process Reference Population

```bash
# QC filtering and LD pruning
plink --bfile [IN_GENO] \
  --keep [reflist] \
  --chr 1-29 \
  --mind 0.1 --maf 0.01 --geno 0.1 \
  --indep-pairwise 50 10 0.2 \
  --make-bed --out REF.1

# Extract pruned SNPs
plink --bfile REF.1 \
  --extract REF.1.prune.in \
  --make-bed --out REF
```

**Output:** `results/[NAME]/Reference.Population/REF.{bed,bim,fam}`

#### Step 2: Create Allele Reference

```bash
# Build allele template (SNP, A1, A2)
awk '{print $2, $5, $6}' REF.bim > alleles.template
```

This ensures mixed population uses same allele coding.

#### Step 3: Process Mixed Population

```bash
# Extract same SNP set with consistent allele coding
plink --bfile [IN_GENO] \
  --keep [mixedlist] \
  --chr 1-29 \
  --mind 0.1 --geno 0.1 \
  --extract REF.1.prune.in \
  --a1-allele alleles.template 3 1 \
  --make-bed --out [NAME]
```

**Output:** `results/[NAME]/Mixed.population/[NAME].{bed,bim,fam}`

#### Step 4: Select Unrelated Individuals

Runs R script `Selection_of_Unrelated_script.R`:

1. **Convert to GDS format:**
   ```r
   snpgdsBED2GDS("REF.bed", "REF.gds")
   ```

2. **Calculate kinship matrix (KING-robust):**
   ```r
   king <- snpgdsIBDKING(gds, type="KING-robust")
   ```

3. **Run PC-AiR:**
   ```r
   pcair_result <- pcair(gds, kinobj=king$kinship, divobj=king$kinship)
   ```
   - Identifies unrelated set for analysis
   - Calculates PCs accounting for relatedness

4. **Save outputs:**
   - `Unrelated_list.txt` - IDs of unrelated samples
   - `PC_scores_PC1-20.csv` - PC scores for all samples
   - Various `.rds` files for reproducibility

**Output:** `results/[NAME]/Reference.Population/Unrelated_list.txt`

#### Step 5: Build Final Datasets

```bash
# Create unrelated reference matched to mixed population
plink --bfile REF \
  --keep Unrelated_list.txt \
  --extract [NAME].snplist \
  --a1-allele alleles.template 3 1 \
  --make-bed --out REF.unrelated
```

#### Step 6: Sanity Checks

Validates:
- SNP count matches between REF.unrelated and mixed
- SNP order is identical (by chr:pos)
- Allele coding is consistent

**Output:** `results/[NAME]/Reference.Population/REF.unrelated.{bed,bim,fam}`

### Stage 2: ADMIXTURE Analysis (`gbc_admixture.sh`)

#### Step 1: Reference ADMIXTURE

```bash
cd results/[NAME]/ADMIXTURE/Supervised_Unrelated/
admixture REF.unrelated.bed 2 -j[CPUS]
```

**Output:** 
- `REF.unrelated.2.Q` - Ancestry proportions
- `REF.unrelated.2.P` - Allele frequencies per population

#### Step 2: Supervised Mode Setup

```bash
# Copy reference allele frequencies
cp REF.unrelated.2.P [NAME].2.P.in
```

The `.P.in` file constrains allele frequencies in supervised analysis.

#### Step 3: Supervised ADMIXTURE on Mixed

```bash
admixture -P [NAME].bed 2 -j[CPUS]
```

**Output:**
- `[NAME].2.Q` - Ancestry proportions for mixed population

#### Step 4: Combine Results with Sample IDs

```r
# combine_Q_FAM.R
fam <- read.table("[NAME].fam")
q <- read.table("[NAME].2.Q")
combined <- cbind(fam[,1:2], q)
write.csv(combined, "[NAME].combined.ancestry.csv")
```

**Output:** `results/[NAME]/ADMIXTURE/Supervised_Unrelated/[NAME].combined.ancestry.csv`

#### Step 5: Generate Visualizations

**Q-Plot (`Q.PLOT.R`):**
```r
# Stacked bar plot of ancestry proportions
# Highlights unrelated individuals used in analysis
```

**PCA Plot (`PCA.plots.R`):**
```r
# Scatter plot of PC1 vs PC2
# Points colored by ancestry proportion
# Shows population structure overlaid with admixture
```

**Outputs:**
- `results/[NAME]/Figures/admixture.png`
- `results/[NAME]/Figures/PCA_Admixture.png`

## Output Files

### Directory Structure

```
results/[NAME]/
├── Reference.Population/
│   ├── REF.{bed,bim,fam}              # Pruned reference
│   ├── REF.unrelated.{bed,bim,fam}    # Unrelated reference
│   ├── Unrelated_list.txt             # Unrelated IDs
│   ├── PC_scores_PC1-20.csv           # PCA results
│   └── alleles.template               # Allele reference
├── Mixed.population/
│   └── [NAME].{bed,bim,fam}           # Mixed population
├── ADMIXTURE/Supervised_Unrelated/
│   ├── [NAME].2.Q                      # Ancestry proportions
│   ├── [NAME].combined.ancestry.csv    # Q + sample IDs
│   ├── REF.unrelated.2.P               # Reference allele frequencies
│   └── REF.unrelated.2.Q               # Reference ancestry (should be pure)
└── Figures/
    ├── admixture.png                   # Q-plot
    └── PCA_Admixture.png               # PCA with ancestry
```

### Ancestry File Format

**`[NAME].combined.ancestry.csv`:**

```csv
FID,IID,Pop1,Pop2
UF,Animal001,0.75,0.25
UF,Animal002,0.42,0.58
```

Columns:
- `FID`, `IID` - Family and individual IDs
- `Pop1` - Proportion from first reference population (e.g., Angus)
- `Pop2` - Proportion from second reference population (e.g., Brahman)

## Interpretation

### Q-Plot

- **X-axis:** Samples (ordered by ancestry proportion)
- **Y-axis:** Ancestry proportion (0-1)
- **Colors:** Different reference populations
- **Markers:** Unrelated samples (used in supervised mode) are highlighted

### PCA Plot

- **Axes:** PC1 vs PC2 (captures most genetic variation)
- **Points:** Individual animals
- **Colors:** Gradient based on ancestry proportion
- **Clustering:** Reference populations should cluster at extremes
- **Admixed animals:** Positioned between reference clusters

### Ancestry Proportions

| Pop1 | Pop2 | Interpretation |
|------|------|----------------|
| ~1.0 | ~0.0 | Pure Pop1 (or nearly pure) |
| ~0.0 | ~1.0 | Pure Pop2 (or nearly pure) |
| 0.5 | 0.5 | F1 cross or balanced admixture |
| 0.75 | 0.25 | 3/4 Pop1, 1/4 Pop2 |

## Quality Control

### Reference Population Quality

```bash
# Check reference Q values (should be ~1.0 for their breed)
head results/[NAME]/ADMIXTURE/Supervised_Unrelated/REF.unrelated.2.Q

# Expected for pure breeds:
# Angus: 1.00 0.00
# Brahman: 0.00 1.00
```

### Unrelated Sample Selection

```bash
# Number of unrelated individuals
wc -l results/[NAME]/Reference.Population/Unrelated_list.txt

# Should be substantial (>50% of total reference)
# If very low (<20%), consider:
# - Reducing kinship threshold in PC-AiR
# - Using different population
```

### SNP Count

```bash
# Check final SNP count
wc -l results/[NAME]/Reference.Population/REF.unrelated.bim

# Typical ranges:
# - Conservative pruning: 80-150K SNPs
# - Moderate pruning: 40-80K SNPs
# - Aggressive pruning: 20-40K SNPs
```

## Parameter Tuning

### QC Thresholds

| Parameter | Conservative | Moderate | Permissive |
|-----------|-------------|----------|------------|
| `mind` | 0.05 | 0.10 | 0.15 |
| `maf` | 0.05 | 0.01 | 0.005 |
| `geno` | 0.05 | 0.10 | 0.15 |
| LD r² | 0.2 | 0.3 | 0.5 |

### LD Pruning Window

Default: `--indep-pairwise 50 10 0.2`

Adjust for different goals:
- **Faster convergence:** `--indep-pairwise 50 10 0.1` (stricter)
- **More SNPs:** `--indep-pairwise 100 10 0.3` (more permissive)

## Common Issues

### Issue 1: ADMIXTURE Won't Converge

**Symptoms:** ADMIXTURE runs >100 iterations without converging

**Solutions:**
- Increase LD pruning stringency (lower r² threshold)
- Remove SNPs with high missingness
- Verify reference populations are truly distinct

### Issue 2: Reference Animals Show Admixture

**Symptoms:** REF.unrelated.2.Q shows values far from 0 or 1

**Solutions:**
- Verify reference animals are truly purebred
- Check for sample swaps or mislabeling
- Remove animals with unexpected ancestry

### Issue 3: Visualization Script Fails

**Symptoms:** R script errors when generating plots

**Solutions:**
```bash
# Check R packages
R -e "library(ggplot2); library(dplyr)"

# Verify input files exist
ls results/[NAME]/ADMIXTURE/Supervised_Unrelated/[NAME].combined.ancestry.csv
ls results/[NAME]/Reference.Population/PC_scores_PC1-20.csv
```

### Issue 4: SNP Mismatch Between Datasets

**Error:** "SNP order mismatch detected"

**Solutions:**
- Re-run pipeline from Stage 1
- Ensure `alleles.template` used correctly
- Check for coordinate system mismatches

## Usage Example

```bash
# 1. Prepare sample ID files
metadata/IDs/MAB.ID 
metadata/IDs/admix.ID 

# 2. Edit breed_composition.sh with paths
bin/pipelines/breed_composition.sh

# 3. Submit pipeline
Run line by line or you can bash or sbatch
sbatch bin/pipelines/breed_composition.sh

# 4. Monitor jobs
squeue -u $USER

# 5. Check results
ls results/GBC_*/Figures/
cat results/GBC_*/ADMIXTURE/Supervised_Unrelated/*.combined.ancestry.csv
```

## Dependencies

### Software
- **PLINK 1.90b3+** - Genotype processing
- **ADMIXTURE 1.3.0+** - Ancestry inference
- **R 3.6+** with packages:
  - `SNPRelate` - GDS format, genetic analysis
  - `GENESIS` - PC-AiR algorithm
  - `gdsfmt` - Genomic data structure
  - `ggplot2` - Plotting
  - `dplyr`, `tidyr` - Data manipulation

### R Package Installation

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "GENESIS", "gdsfmt"))

# CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr"))
```

## SLURM Resources

Default configuration:

```bash
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=48:00:00
```

Scale based on dataset:
- Small (<5K samples, <100K SNPs): 3 CPUs, 8GB/CPU
- Medium (5-20K samples): 6 CPUs, 12GB/CPU
- Large (>20K samples): 12 CPUs, 16GB/CPU

**Runtime estimates:**
- Population prep: 1-4 hours
- ADMIXTURE: 2-8 hours
- Visualization: <10 minutes

## Next Steps

After breed composition analysis:
- Use ancestry proportions as covariates in GWAS
- Identify purebred vs crossbred animals
- Correct for population stratification
- Design breeding strategies
- Validate pedigree information

## References

Relevant methods:
- **ADMIXTURE:** Alexander et al. (2009) Genome Research
- **PC-AiR:** Conomos et al. (2015) Genetic Epidemiology
- **KING-robust:** Manichaikul et al. (2010) Bioinformatics
