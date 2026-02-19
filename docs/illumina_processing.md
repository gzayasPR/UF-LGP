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
ALLELE_POS=4                                        # Column for alleles (4=Forward, 8=AB)
```

## Input Files

### 1. Illumina Metadata CSV (`Illumina_meta_csv`)

Lists all FinalReport files to process.

**Format:** CSV with header `full_path`
```csv
full_path
/blue/mateescu/raluca/Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt
/blue/mateescu/raluca/Univ_of_Florida_Mateescu_BOVF250V1_20160401_FinalReport.txt
```

- Each row = one FinalReport file
- Full absolute paths required
- No blank lines

### 2. FinalReport Files

Illumina's standard genotyping output format containing raw genotype calls.

**File Structure:**
```
[Header]
GSGT Version    2.0.4
Processing Date 11/9/2021 9:43 AM
Content         GGP-F250_20000778_A1.bpm
Num SNPs        221115
Total SNPs      221115
Num Samples     432
Total Samples   432
[Data]
SNP Name        Sample ID       Allele1 - Forward       Allele2 - Forward       Allele1 - Top   Allele2 - Top   Allele1 - AB    Allele2 - AB    GC Score        X       Y
1-1001021-A-C-rs439448877       Z000812288      A       A       A       A       A       A       0.4616  0.675   0.020
1-100158803-T-G Z000812288      T       T       A       A       A       A       0.2153  0.112   0.047
1-100158806-C-G Z000812288      C       C       G       G       B       B       0.1872  0.104   0.564
```

**Key Sections:**

- **`[Header]`**: Metadata about the genotyping run
  - `Content`: SNP chip version (e.g., GGP-F250 = GeneSeek Genomic Profiler 250K)
  - `Num SNPs`: Number of markers on the chip
  - `Num Samples`: Animals genotyped in this run
  
- **`[Data]`**: Genotype calls (one row per SNP per sample)

**Important Columns:**

| Column | Description | Example | Used By Pipeline |
|--------|-------------|---------|------------------|
| `SNP Name` | Marker identifier | `1-1001021-A-C-rs439448877` | YES - SNP ID |
| `Sample ID` | Animal identifier | `Z000812288` | YES - Individual ID |
| `Allele1 - Forward` | First allele (forward strand) | `A` | YES (if ALLELE_POS=4) |
| `Allele2 - Forward` | Second allele (forward strand) | `C` | YES (if ALLELE_POS=4) |
| `Allele1 - Top` | First allele (TOP strand) | `A` | Optional |
| `Allele2 - Top` | Second allele (TOP strand) | `C` | Optional |
| `Allele1 - AB` | First allele (A/B coding) | `A` | YES (if ALLELE_POS=8) |
| `Allele2 - AB` | Second allele (A/B coding) | `B` | YES (if ALLELE_POS=8) |
| `GC Score` | Genotype quality score (0-1) | `0.4616` | Used for QC |
| `X` | Normalized intensity (allele 1) | `0.675` | Used for QC |
| `Y` | Normalized intensity (allele 2) | `0.020` | Used for QC |

**Allele Coding Systems:**

- **Forward Strand**: Actual DNA sequence orientation (A, C, G, T)
- **TOP Strand**: Illumina's standardized orientation (A, C, G, T)
- **AB Coding**: Arbitrary labels (A or B) - chip-specific, not transferable

> **Note:** This pipeline uses **Forward strand alleles** (`ALLELE_POS=4`) by default to maintain consistency with reference genome coordinates. AB coding should only be used if you need chip-specific analysis.

### 3. SNP Map File (`SNP_MAP`)

Maps Illumina SNP names to chromosome and position.

**Format:** Tab-delimited, no header
```
SNP_Name        Chromosome      Position
1-1001021-A-C-rs439448877       1       1001021
1-100158803-T-G 1       100158803
1-100158806-C-G 1       100158806
```

- Column 1: Illumina SNP identifier (matches FinalReport `SNP Name`)
- Column 2: Chromosome number (1-29, X)
- Column 3: Base pair position

### 4. ID Mapping File (`UPDATE_IDS`)

Standardizes animal IDs from Illumina lab format to project-specific identifiers.

**Format:** Tab-delimited with header
```
DateRec  IlluminaID_No_Spaces    PID_No_Spaces     UFID_No_Spaces
20160331 V_2130309_Plate5_A01    MAB_BRA           2130309
20160401 3060076                 MAB_BRA           3060076
20171127 950514                  Brahman           950514
20220505 Dc_529                  Gonella           Dc_529
20171228 500163                  SeminThermo       500163
```

**Column Definitions:**

- **`DateRec`** (OldFID): Processing date from Illumina FinalReport header
  - Format: YYYYMMDD
  - Found in FinalReport file for example:
      - Univ_of_Florida_Marsella_BOVF250V1_20211108_FinalReport.txt (Date:20211108)
      - Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt (Date:20160331)
  - Used as temporary Family ID during conversion
  
- **`IlluminaID_No_Spaces`** (OldIID): Sample ID from Illumina FinalReport
  - Exact match to `Sample ID` column in FinalReport `[Data]` section
  - May include plate positions (e.g., `V_2130309_Plate5_A01`) or plain IDs
  - Illumina IDs may have multiple spaces like V 2130309 Plate5 A01 or Dc 529. 
      - These spaces will be replaced in script with "_". 
      - Thus when renaming no spaces are permitted in the ID mapping file. 
  - Lab-assigned identifiers, not standardized across projects
  
- **`PID_No_Spaces`** (NewFID): Target project/breed identifier
  - Examples: `MAB_BRA`, `Brahman`, `Gonella`, `SeminThermo`
  - Groups animals by breeding program or research project
  - Used as Family ID in final PLINK files
  
- **`UFID_No_Spaces`** (NewIID): University of Florida animal ID
  - Unique animal identifier
  - Standardized across all genotyping runs
  - Used as Individual ID in final PLINK files

**Mapping Logic:**
```
Illumina FinalReport                    →    Final PLINK Format
─────────────────────────────────────────    ──────────────────
Processing Date: 03/31/2016                  FID: MAB_BRA
Sample ID: V_2130309_Plate5_A01        →    IID: 2130309
```

**Why This Mapping is Needed:**

1. **Date-based batching**: Illumina processes samples in batches by date, not by project
2. **Lab vs. farm IDs**: Lab uses plate positions; farm uses permanent animal IDs  
3. **Cross-run consistency**: Same animal may be genotyped multiple times with different lab IDs
4. **Project organization**: Groups animals by biological/breeding relevance, not processing date

> **Important:** Every unique combination of `DateRec` + `IlluminaID_No_Spaces` in your FinalReports should have a corresponding entry in this file (even if ID is the same), if not you can find them easily in your plink.fam file. 

### 5. Chromosome Update File (`CHR_UPD`)

Corrects chromosome assignments to match reference genome.

**Format:** Tab-delimited, no header
```
SNP_Name        Chromosome
1-1001021-A-C-rs439448877       1
1-100158803-T-G 1
```

- Used when original SNP map has outdated or incorrect chromosome assignments
- Ensures compatibility with ARS-UCD1.2 reference

### 6. Base Pair Update File (`BP_UPD`)

Updates SNP positions to ARS-UCD1.2 reference genome coordinates.

**Format:** Tab-delimited, no header
```
SNP_Name        Position
1-1001021-A-C-rs439448877       1001021
1-100158803-T-G 100158803
```

- Critical for accurate variant mapping
- Positions must match current reference genome build

### 7. Indel Reference File (`INDELS`)

Specifies correct allele coding for insertion/deletion variants.

**Format:** Tab-delimited
```
SNP_Name        RefAllele       AltAllele
rs123456        A       -
rs789012        -       CGTA
```

- `-` indicates deletion
- Multi-base strings indicate insertions
- Ensures proper variant representation

## Quality Control Notes

**FinalReport Quality Indicators:**

- **GC Score** < 0.15: Low confidence genotype (may be excluded)
- **X ≈ 0 and Y ≈ 0**: No-call/missing genotype
- **High X, low Y**: Homozygous for allele 1
- **Low X, high Y**: Homozygous for allele 2
- **Moderate X and Y**: Heterozygous

**Common Issues:**

1. **Mismatched SNP counts**: Check if FinalReport contains all expected markers
2. **Sample ID format changes**: Verify `UPDATE_IDS` file is current
3. **Coordinate mismatches**: Ensure `BP_UPD` and `CHR_UPD` match reference genome version
4. **Allele flips**: Pipeline handles strand flips automatically using Forward strand

---

> **Important:** Always verify your `ALLELE_POS` setting matches your analysis needs. Forward strand (column 4) is recommended for compatibility with reference genomes. If using only this you can use Illumina alleles (A/B) which is column 8. 

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
