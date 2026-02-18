# Illumina Array Processing Pipeline

Convert Illumina Bovine 250K FinalReport files to PLINK format with proper genome coordinates.

## Overview

This pipeline processes multiple Illumina FinalReport files in batch, performing:
1. File validation and format checking
2. Genotype extraction and ID standardization
3. LGEN → PLINK binary conversion
4. Coordinate mapping to ARS-UCD1.2 reference genome
5. Quality control and error detection

## Pipeline Entry Point

**Script:** `bin/pipelines/illumina_2_plink.sh`

## Required Configuration

Edit the following variables in `illumina_2_plink.sh`:

```bash
proj_env="/path/to/bin/project_env.sh"
NAME="250K_2026"                                    # Output name
Illumina_meta_csv="/path/to/metadata/finalreports.csv"
SNP_MAP="/path/to/metadata/SNP_Map.txt"
UPDATE_IDS="/path/to/metadata/Feb_2024.txt"        # ID mapping file
CHR_UPD="/path/to/metadata/ARS.UF.Chr.Update.txt"  # Chromosome updates
BP_UPD="/path/to/metadata/ARS.BP.2024.txt"         # Base pair positions
INDELS="/path/to/metadata/INDELS.ID"               # Allele reference
ALLELE_POS=4                                        # Column for alleles (4=real, 8=Illumina)
```

## Input Files

### 1. Illumina Metadata CSV (`Illumina_meta_csv`)

**Format:** CSV with header `full_path`

```csv
full_path
/blue/mateescu/raluca/Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt
/blue/mateescu/raluca/Univ_of_Florida_Mateescu_BOVF250V1_20160401_FinalReport.txt
```

### 2. FinalReport Files

Standard Illumina format with 10-line header:
- Line 1-9: Metadata (chip info, sample count, etc.)
- Line 10: Column headers
- Line 11+: Genotype data (11 columns expected)

**Expected Columns:**
```
SNP_Name  Sample_ID  Allele1_Top  Allele2_Top  GC_Score  X  Y  Chr  Position  ...
```

### 3. Reference Files

**SNP_Map.txt:** Maps SNP names to ARS-UCD1.2 coordinates
```
SNP_Name    Chr    Position    Allele1    Allele2
ARS-BFGL-BAC-10172    1    123456    A    G
```

**UPDATE_IDS:** Maps Sample_ID to standardized IDs
```
OldFID  OldIID  NewFID  NewIID
20160331    Sample001    UFID    Animal001
```

**CHR_UPD:** Chromosome name corrections
```
SNP_Name    Chr
ARS-BFGL-BAC-10172    1
```

**BP_UPD:** Base pair position updates
```
SNP_Name    Position
ARS-BFGL-BAC-10172    123456
```

**INDELS.ID:** Allele corrections (INDEL normalization)
```
SNP_Name    Allele1    Allele2
Hapmap43437-BTA-101873    A    G
```

## Processing Steps

### Stage 1: File Preparation (`illumina_prep.sh`)

1. **Copy files** from source directories to working directory
2. **Validate format:**
   - Check header row count (should be 10)
   - Verify column count (should be 11)
   - Check "Num Samples" field
3. **Detect format errors** (e.g., extra spaces in Sample_ID column)

### Stage 2: PLINK Conversion (`make_ill_plink.sh`)

1. **Combine FinalReports:**
   - Remove 10-line headers from each file
   - Concatenate into `all.txt`

2. **Extract genotypes:**
   - Parse alleles from column `ALLELE_POS` (4 or 8)
   - Build LGEN format (FID, IID, SNP, Allele1, Allele2)
   - Replace spaces in IDs with underscores
   - Convert missing data (`- -`) to PLINK format (`0 0`)

3. **Create PLINK files:**
   ```bash
   plink --lgen [NAME].lgen --fam [NAME].fam --map [NAME].map --make-bed
   ```

4. **Update IDs:**
   ```bash
   plink --bfile [NAME] --update-ids [UPDATE_IDS] --make-bed
   ```

5. **Update genomic coordinates:**
   ```bash
   plink --bfile UFID_[NAME] --update-chr [CHR_UPD] --make-bed
   plink --bfile mydata2 --update-map [BP_UPD] --make-bed
   plink --bfile mydata3 --update-alleles [INDELS] --make-bed
   ```

## Output Files

### Final PLINK Binary Files

**Location:** `results/[NAME]/`

```
UFID_[NAME]_Illumina.ID_ARS.bed    # Binary genotype data
UFID_[NAME]_Illumina.ID_ARS.bim    # SNP information (6 columns)
UFID_[NAME]_Illumina.ID_ARS.fam    # Sample information (6 columns)
```

### Intermediate Files

**Location:** `results/[NAME]/Make_PLINK/`

- `all.txt` - Combined FinalReports (no headers)
- `Geno.all.txt` - Extracted genotypes
- `ID.txt` - Date + Sample_ID mapping
- `[NAME].lgen` - PLINK long format
- `[NAME].map` - SNP map
- `[NAME].fam` - Sample pedigree

### Copied FinalReports

**Location:** `results/[NAME]/Illumina/`

Original FinalReport files copied for processing.

## SLURM Configuration

Default resource allocation:

```bash
#SBATCH --account=mateescu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=12gb
```

Adjust based on dataset size:
- Small (<1000 samples): 8 CPUs, 8GB/CPU
- Medium (1000-5000): 10 CPUs, 12GB/CPU
- Large (>5000): 16 CPUs, 16GB/CPU

## Quality Control

### Automatic Checks

1. **File existence:** All input files must exist
2. **Header validation:** 10 header lines in FinalReports
3. **Column count:** 11 columns in genotype section
4. **Sample count:** Match between "Num Samples" and actual rows
5. **PLINK validation:** Automatic checks during conversion

### Manual Verification

After pipeline completion:

```bash
# Check output dimensions
plink --bfile results/[NAME]/UFID_[NAME]_Illumina.ID_ARS --freq

# Verify coordinate mapping
head results/[NAME]/UFID_[NAME]_Illumina.ID_ARS.bim

# Check missing data
plink --bfile results/[NAME]/UFID_[NAME]_Illumina.ID_ARS --missing
```

## Common Issues

### Issue 1: Column Count Mismatch

**Error:** `head -20 ${file} |awk '{print NF}' returns 12 instead of 11`

**Cause:** Extra space in Sample_ID column (e.g., "Dc 14" instead of "Dc_14")

**Solution:** Use sed to replace spaces:
```bash
sed -i 's/Dc 14/Dc_14/g' [FinalReport.txt]
```

### Issue 2: ALLELE_POS Setting

**Error:** Unexpected allele values in output

**Solution:**
- Use `ALLELE_POS=4` for real genotypes (Allele1_Top, Allele2_Top)
- Use `ALLELE_POS=8` for Illumina-specific format

### Issue 3: Missing Reference Files

**Error:** Cannot find SNP_Map.txt or UPDATE_IDS

**Solution:** Check paths in `metadata/Illumina_to_Plink/` directory

## Usage Example

```bash
# 1. Prepare metadata CSV
cat > metadata/250K_finalreports.csv << EOF
full_path
/blue/genotypes/Univ_of_Florida_20160331_FinalReport.txt
/blue/genotypes/Univ_of_Florida_20160401_FinalReport.txt
EOF

# 2. Edit pipeline script with correct paths
vim bin/pipelines/illumina_2_plink.sh

# 3. Submit job
sbatch bin/pipelines/illumina_2_plink.sh

# 4. Monitor progress
tail -f bash_out/Illumina_2_PLINK/makePLINK_250K_2026_*.log

# 5. Check outputs
ls -lh results/250K_2026/UFID_250K_2026_Illumina.ID_ARS.*
```

## Dependencies

- **PLINK 1.90b3+** - Format conversion and QC
- **GNU coreutils** - awk, sed, paste, sort
- **dos2unix** - Line ending conversion (optional)

## Performance Notes

- **Runtime:** ~1-4 hours depending on sample count
- **Memory:** Scales with SNP count (~200-250K SNPs = 8-12GB)
- **Disk:** Intermediate files can be large (5-20GB total)

## Next Steps

After successful processing, the PLINK files can be used for:
- Quality control filtering (`--mind`, `--geno`, `--maf`)
- Merging with other datasets
- Population structure analysis
- GWAS or genomic prediction

See other pipeline documentation for downstream analysis options.
