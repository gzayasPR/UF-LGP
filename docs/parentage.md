# Parentage Verification Pipeline

KING-based relatedness inference for identifying parent-offspring relationships.

## Overview

This pipeline uses KING (Kinship-based INference for Gwas) to identify parent-offspring pairs through genome-wide relatedness estimation. It processes a list of query individuals and identifies their parents from a reference population.

## Pipeline Entry Point

**Script:** `bin/pipelines/check_parentage.sh`

## Required Configuration

Edit these variables in `check_parentage.sh`:

```bash
proj_env="/path/to/bin/project_env.sh"
NAME="Parentage_2026"
IN_GENO="/path/to/genotypes"        # PLINK fileset
reflist="/path/to/MAB.ID"           # Reference population (potential parents)
INDVI_LIST="/path/to/check.parentage.ID"  # Individuals to verify

# QC parameters
mind=0.1    # Max individual missingness
maf=0.01    # Minor allele frequency  
geno=0.1    # Max SNP missingness
```

## Input Files

### 1. Genotype Data (`IN_GENO`)

**Format:** PLINK binary fileset (`.bed/.bim/.fam`)

**Requirements:**
- Must include both query individuals AND potential parents
- Autosomes 1-29 recommended (sex chromosomes can be included)
- Good quality genotypes (high call rate)

### 2. Reference Population (`reflist`)

**Format:** Two-column space/tab-delimited (FID, IID)

```
Family1  PotentialDad001
Family1  PotentialMom001
Family2  PotentialDad002
```

**Purpose:** Pool of potential parents to search within

**Best practices:**
- Include all possible parents
- Include siblings if uncertain
- Larger reference = better chance of finding parents

### 3. Query Individuals (`INDVI_LIST`)

**Format:** One ID per line (IID only, no FID)

```
2190762
2210676
2210764
```

**Requirements:**
- Must exist in the genotype data
- Can overlap with reference population

## Processing Steps

### Step 1: Data Preparation and Pruning

```bash
# Filter and prune data
plink --bfile [IN_GENO] \
  --keep [reflist] \
  --chr 1-29 \
  --maf 0.05 --geno 0.01 --mind 0.1 \
  --hwe 0.00001 \
  --indep-pairwise 5000 10 0.5 \
  --make-bed --out Temp

# Extract pruned SNPs
plink --bfile Temp \
  --extract Temp.prune.in \
  --make-bed --out king_input
```

**Rationale:**
- **Large window (5000):** Captures long-range LD patterns
- **High r² threshold (0.5):** Retains more SNPs for accurate kinship
- **HWE filter:** Removes genotyping errors
- **Strict missingness:** Ensures accurate IBS0 calculation

### Step 2: KING Relatedness Analysis

```bash
king -b king_input.bed --related --sexchr 30
```

**Outputs:**
- `king.kin` - Pairwise relatedness for all samples
- `king.kin0` - Sparse format (closely related pairs only)

### Step 3: Extract Per-Individual Results

For each individual in `INDVI_LIST`:

```bash
# Extract all relationships for this individual
grep [INDVI] king.kin > parentage_results/[INDVI].kin

# Add header
sed -i "1s/^/[HEADER]\n/" parentage_results/[INDVI].kin
```

### Step 4: Identify Parent-Offspring Pairs

```bash
# Filter for InfType = "PO" (Parent-Offspring)
awk 'NR>1 && $15=="PO"' parentage_results/[INDVI].kin \
  > parentage_results/parent_offspring/[INDVI]_parents.txt
```

### Step 5: Generate Summary Report

```bash
# Create CSV summary
echo "Individual,Parent,InfType,Kinship,IBS0" > parent_offspring_summary.csv

# Extract parent IDs and statistics
for each individual:
  awk '{print indv, parent_id, $15, $11, $8}' >> parent_offspring_summary.csv
```

## Output Files

### Directory Structure

```
results/Parentage_2026/
├── parentage_results/
│   ├── 2190762.kin                    # All relationships for individual
│   ├── 2210676.kin
│   ├── all_individuals.kin            # Combined results
│   └── parent_offspring/
│       ├── 2190762_parents.txt        # Only PO relationships
│       └── 2210676_parents.txt
├── parent_offspring_summary.csv       # CSV summary
└── Temp.irem                          # Excluded individuals (high missingness)
```

### Parent-Offspring Summary CSV

**Format:**
```csv
Individual,Parent,InfType,Kinship,IBS0
2190762,1040084,PO,0.2456,0.0012
2210676,2000348,PO,0.2501,0.0008
```

**Columns:**
- `Individual` - Query animal ID
- `Parent` - Identified parent ID
- `InfType` - Relationship type ("PO" for parent-offspring)
- `Kinship` - Kinship coefficient (~0.25 for parent-offspring)
- `IBS0` - Proportion of loci with zero alleles shared IBS (~0 for true parent-offspring)

## KING Output Interpretation

### Key Columns in .kin Files

| Column | Field | Description | PO Value |
|--------|-------|-------------|----------|
| 2 | ID1 | First individual | Query or Parent |
| 3 | ID2 | Second individual | Parent or Query |
| 8 | IBS0 | Proportion SNPs with 0 shared alleles | ≈0 |
| 11 | Kinship | Kinship coefficient | ≈0.25 |
| 15 | InfType | Inferred relationship | "PO" |

### Relationship Types (InfType)

| InfType | Relationship | Kinship | IBS0 | Notes |
|---------|-------------|---------|------|-------|
| **PO** | Parent-Offspring | ~0.25 | ~0 | Target relationship |
| FS | Full-Sibling | ~0.25 | >0 | Same parents, IBS0 > 0 |
| 2nd | Second-degree | ~0.125 | - | Half-sibs, grandparent, etc. |
| 3rd | Third-degree | ~0.0625 | - | Cousins |
| UN | Unrelated | <0.05 | - | No close relationship |
| Dup/MZ | Duplicate/Twins | ~0.5 | ~0 | Same individual or MZ twins |

### Distinguishing Parent-Offspring from Full-Siblings

**Parent-Offspring:**
- Kinship ≈ 0.25
- **IBS0 ≈ 0** (parents share at least one allele at every locus)
- InfType = "PO"

**Full-Siblings:**
- Kinship ≈ 0.25
- **IBS0 > 0** (siblings can have 0 shared alleles at some loci)
- InfType = "FS"

**Key discriminator:** IBS0 is near-zero for true parent-offspring pairs.

## Quality Control

### Genotype Quality

```bash
# Check individual missingness
plink --bfile king_input --missing

# Individuals with >10% missing data are excluded
cat Temp.irem
```

### Relationship Statistics

**Expected values for true parent-offspring:**
- Kinship: 0.23 - 0.27 (theoretical 0.25)
- IBS0: 0.000 - 0.005 (close to 0)

**Red flags:**
- Kinship < 0.20 or > 0.30: Possible genotyping error
- IBS0 > 0.01: May not be true parent-offspring
- No "PO" relationships found: Parents not in reference panel

### Verification Steps

```bash
# Count parent-offspring pairs found
grep "PO" parent_offspring_summary.csv | wc -l

# Check kinship distribution for PO pairs
awk -F',' 'NR>1 {print $4}' parent_offspring_summary.csv | \
  awk '{sum+=$1; sumsq+=$1*$1; n++} END {print "Mean:", sum/n, "SD:", sqrt(sumsq/n - (sum/n)^2)}'

# Should be: Mean ≈ 0.25, SD < 0.02
```

## Parameter Tuning

### LD Pruning

Default: `--indep-pairwise 5000 10 0.5`

**Adjust for:**
- **More SNPs (slower, more accurate):** `--indep-pairwise 10000 10 0.6`
- **Fewer SNPs (faster):** `--indep-pairwise 1000 10 0.3`

### QC Thresholds

| Parameter | Stringent | Default | Permissive |
|-----------|-----------|---------|------------|
| `mind` | 0.05 | 0.10 | 0.15 |
| `maf` | 0.05 | 0.05 | 0.01 |
| `geno` | 0.005 | 0.01 | 0.05 |
| HWE | 1e-6 | 1e-5 | 1e-4 |

## Common Scenarios

### Scenario 1: Both Parents in Reference

**Result:** Two PO relationships identified

```csv
Individual,Parent,InfType,Kinship,IBS0
2190762,1040084,PO,0.2456,0.0012
2190762,1040085,PO,0.2489,0.0008
```

### Scenario 2: One Parent in Reference

**Result:** One PO relationship identified

```csv
Individual,Parent,InfType,Kinship,IBS0
2210676,2000348,PO,0.2501,0.0008
```

### Scenario 3: No Parents in Reference

**Result:** No PO relationships, may see other relationships

```
✗ No parent-offspring relationships found for 2210764
```

**Possible explanations:**
- Parents not genotyped
- Parents not included in `reflist`
- Sample swap or mislabeling

### Scenario 4: Unexpected Relationships

**Multiple PO pairs (>2):**
- Check for duplicates or identical twins in reference
- Verify sample IDs

**FS instead of PO:**
- Individual is a full-sibling, not offspring
- Update pedigree accordingly

## Common Issues

### Issue 1: No Relationships Found

**Symptoms:** Empty parent_offspring results for all individuals

**Solutions:**
1. Verify query IDs exist in genotype data:
   ```bash
   awk '{print $2}' king_input.fam | grep -f [INDVI_LIST]
   ```

2. Check reference population size:
   ```bash
   wc -l [reflist]
   # Should have potential parents
   ```

3. Verify data quality:
   ```bash
   plink --bfile king_input --missing
   ```

### Issue 2: Kinship Values Outside Expected Range

**Symptoms:** Kinship far from 0.25 for "PO" relationships

**Possible causes:**
- Genotyping errors
- Sample contamination
- Platform differences (array vs sequencing)

**Solutions:**
- Increase QC stringency
- Check for batch effects
- Verify sample identities

### Issue 3: High IBS0 for "PO" Pair

**Symptoms:** IBS0 > 0.01 for InfType="PO"

**Interpretation:**
- May not be true parent-offspring
- Could be full-siblings misclassified
- Genotyping error

**Action:**
- Manually review relationship
- Check pedigree records
- Consider other evidence (birth dates, etc.)

### Issue 4: Memory Errors

**Symptoms:** KING crashes or "Killed" message

**Solutions:**
- Reduce SNP count with stricter LD pruning
- Process in batches
- Increase memory allocation:
  ```bash
  #SBATCH --mem-per-cpu=16gb
  ```

## Usage Example

```bash
# 1. Create list of individuals to verify
cat > metadata/check.parentage.ID << EOF
2190762
2210676
2210764
EOF

# 2. Verify reference population list exists
head metadata/IDs/MAB.ID

# 3. Edit check_parentage.sh with correct paths
vim bin/pipelines/check_parentage.sh

# 4. Submit job
sbatch bin/pipelines/check_parentage.sh

# 5. Monitor progress
tail -f bash_out/Parentage_2026/KING_*.log

# 6. Check results
cat results/Parentage_2026/parent_offspring_summary.csv

# 7. Review specific individual
cat results/Parentage_2026/parentage_results/parent_offspring/2190762_parents.txt
```

## Batch Processing

For large numbers of individuals:

```bash
# Split into batches of 100
split -l 100 all_individuals.list batch_

# Submit multiple jobs
for batch in batch_*; do
  sbatch --export=INDVI_LIST=$batch bin/pipelines/check_parentage.sh
done

# Combine results
cat results/Parentage_*/parent_offspring_summary.csv > all_parentage_results.csv
```

## Dependencies

- **PLINK 1.90b3+** - Genotype processing
- **KING 2.3.0+** - Relatedness inference
- **Standard Unix tools** - grep, awk, sed

## SLURM Resources

Default configuration:

```bash
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=48:00:00
```

**Runtime estimates:**
- Small (<5K samples, <100K SNPs): 1-2 hours
- Medium (5-20K samples): 3-6 hours
- Large (>20K samples): 6-24 hours

**Scaling:**
- Memory scales with samples²
- Runtime scales with samples × SNPs

## Validation

### Cross-Reference with Pedigree

```bash
# If pedigree data available, compare results
# Create pedigree file: offspring, sire, dam
join -1 1 -2 1 <(sort pedigree.txt) <(sort parent_offspring_summary.csv)
```

### Check for Mendelian Inconsistencies

```bash
# After verifying parentage, test for inheritance errors
plink --bfile [genotypes] \
  --mendel \
  --set-me-missing \
  --mendel-duos \
  --make-bed --out verified_parentage
```

## Next Steps

After parentage verification:
- Update pedigree records
- Identify pedigree errors or sample swaps
- Use verified relationships for:
  - Heritability estimation
  - Genomic prediction (relationship matrix)
  - Imputation (family-based)
  - Quality control (remove duplicates/errors)

## References

- **KING:** Manichaikul et al. (2010) "Robust relationship inference in genome-wide association studies." Bioinformatics 26(22):2867-2873.
- **IBS0:** Anderson et al. (2010) "Data quality control in genetic case-control association studies." Nature Protocols 5:1564-1573.
