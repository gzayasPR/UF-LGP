# File Formats Reference

Documentation of input and output file formats used in UF-LGP pipelines.

## Input File Formats

### Sample ID Lists

**Usage:** Specifying individuals for analysis

**Format:** Two-column space or tab-delimited
```
FID    IID
Family1    Animal001
Family1    Animal002
UF    2190762
```

**Columns:**
- `FID` - Family ID (can be arbitrary, often set to project name)
- `IID` - Individual ID (must be unique within FID)

**Notes:**
- No header row
- Space or tab delimited (not comma)
- No quotes unless IDs contain spaces (not recommended)

---

### PLINK Binary Format (.bed/.bim/.fam)

**Usage:** Genotype storage and analysis

#### .bed File (Binary Genotype Data)
- Binary file containing packed genotypes
- 2 bits per genotype: 00=homozygous ref, 01=missing, 10=heterozygous, 11=homozygous alt
- Not human-readable

#### .bim File (SNP Information)
**Format:** 6 columns, tab-delimited, no header
```
1    ARS-BFGL-BAC-10172    0    123456    A    G
1    ARS-BFGL-BAC-10176    0    156789    C    T
```

**Columns:**
1. Chromosome code
2. Variant ID
3. Position in morgans (usually 0)
4. Base-pair position
5. Allele 1 (reference)
6. Allele 2 (alternate)

#### .fam File (Sample Information)
**Format:** 6 columns, space-delimited, no header
```
Family1    Animal001    0    0    0    -9
Family1    Animal002    0    0    0    -9
```

**Columns:**
1. Family ID
2. Individual ID
3. Paternal ID (0 = missing)
4. Maternal ID (0 = missing)
5. Sex (1=male, 2=female, 0=unknown)
6. Phenotype (-9=missing, 1=control, 2=case, or quantitative value)

---

### VCF Format (.vcf / .vcf.gz)

**Usage:** Variant call storage with sample genotypes

**Header lines:** Start with `##`
```
##fileformat=VCFv4.2
##contig=<ID=1,length=158534110>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">
```

**Column header:** Starts with `#CHROM`
```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
```

**Data lines:**
```
1  123456  rs123  A  G  .  PASS  .  GT:GP  0/1:0.05,0.90,0.05  1/1:0.01,0.09,0.90
```

**KEY FIELDS:**
- `GT` - Genotype (0/0, 0/1, 1/1, ./.)
- `GP` - Genotype probabilities (used in Skim-Seek filtering)
- `FILTER` - PASS or LOWCONF (used in Skim-Seek)

**Compression:**
- `.vcf.gz` must be BGZF-compressed (use bcftools view -Oz)
- Requires tabix index (`.tbi`) for random access

---

### Illumina FinalReport Format

**Usage:** Illumina SNP array raw output

**Structure:**
- Lines 1-9: Metadata header
- Line 10: Column names
- Lines 11+: Genotype data

**Header example:**
```
[Header]
GSGT Version    2.0.4
Processing Date    3/31/2016 9:38 AM
Content        BovineSNP50_v3_A1.bpm
Num SNPs    54609
Num Samples    96
[Data]
SNP Name    Sample ID    Allele1 - Top    Allele2 - Top    GC Score    X    Y    Chr    Position    GT Score    ...
```

**Data columns (11 total):**
1. SNP Name
2. Sample ID
3. Allele1 - Top
4. Allele2 - Top
5. GC Score (genotype quality)
6. X (intensity)
7. Y (intensity)
8. Chr
9. Position
10. GT Score
11. (varies)

**Notes:**
- Real genotypes in columns 3-4 (ALLELE_POS=4)
- Alternative format uses columns 7-8 (ALLELE_POS=8)
- Missing data: `- -` or empty

---

## Metadata File Formats

### Illumina Metadata CSV

**Usage:** List of FinalReport files to process

**Format:** CSV with header
```csv
full_path
/blue/data/Univ_of_Florida_BOVF250V1_20160331_FinalReport.txt
/blue/data/Univ_of_Florida_BOVF250V1_20160401_FinalReport.txt
```

**Requirements:**
- Must have `full_path` column header
- One file per row
- Full absolute paths

---

### Skim-Seek Metadata CSV

**Usage:** Sample information for VCF merging

**Format:** CSV with 5 required columns
```csv
path,SkimSeek_ID,ProjectID,IID,NewID
/blue/data/sample1.vcf.gz,SSK001,Project1,Animal001,UF_Animal001
/blue/data/sample2.vcf.gz,SSK002,Project1,Animal002,UF_Animal002
```

**Columns:**
- `path` - Full path to per-sample VCF
- `SkimSeek_ID` - Sample name in VCF header (must match exactly)
- `ProjectID` - Project identifier (optional, for tracking)
- `IID` - Original individual ID (optional)
- `NewID` - Desired final sample name in merged VCF

**Requirements:**
- Header row required
- No duplicate `SkimSeek_ID` or `NewID` values
- All VCF files must exist at specified paths

---

### ID Update Files

**Usage:** Renaming samples or updating coordinates

#### PLINK --update-ids Format
**Format:** 4 columns, space-delimited, no header
```
OldFID    OldIID    NewFID    NewIID
20160331    Sample001    UFID    Animal001
20160331    Sample002    UFID    Animal002
```

#### PLINK --update-chr Format
**Format:** 2 columns, space-delimited, no header
```
SNP_Name    Chr
ARS-BFGL-BAC-10172    1
ARS-BFGL-BAC-10176    1
```

#### PLINK --update-map Format
**Format:** 2 columns, space-delimited, no header
```
SNP_Name    Position
ARS-BFGL-BAC-10172    123456
ARS-BFGL-BAC-10176    156789
```

#### PLINK --update-alleles Format
**Format:** 3 columns, space-delimited, no header
```
SNP_Name    Allele1    Allele2
ARS-BFGL-BAC-10172    A    G
ARS-BFGL-BAC-10176    C    T
```

---

## Output File Formats

### ADMIXTURE Q File (.Q)

**Usage:** Ancestry proportions per individual

**Format:** Space-delimited, no header, one row per sample
```
0.9876    0.0124
0.0023    0.9977
0.4567    0.5433
```

**Interpretation:**
- Columns represent ancestral populations (K populations)
- Rows correspond to samples (same order as .fam file)
- Values sum to 1.0 per row
- Example above is K=2 (two populations)

---

### ADMIXTURE P File (.P)

**Usage:** Allele frequencies per population

**Format:** Space-delimited, one row per SNP
```
0.123    0.876
0.456    0.544
```

**Interpretation:**
- Columns represent ancestral populations
- Rows correspond to SNPs (same order as .bim file)
- Frequency of A1 allele in each population

---

### Combined Ancestry CSV

**Usage:** Q matrix with sample IDs

**Format:** CSV with header
```csv
FID,IID,Pop1,Pop2
UF,Animal001,0.9876,0.0124
UF,Animal002,0.4567,0.5433
```

**Columns:**
- `FID`, `IID` - Sample identifiers
- `Pop1`, `Pop2` - Ancestry proportions (one column per K)

---

### KING .kin File

**Usage:** Pairwise relatedness estimates

**Format:** Space-delimited with header
```
FID    ID1    ID2    N_SNP    Z0    Phi    HetHet    IBS0    HetConc    HomIBS0    Kinship    IBD1Seg    IBD2Seg    PropIBD    InfType    Error
FAM    A001    A002    50000    0.25    0.25    0.15    0.001    0.90    0.001    0.2456    5    2    0.75    PO    0
```

**Key columns:**
- `ID1`, `ID2` - Sample pair
- `Kinship` - Kinship coefficient
- `IBS0` - Proportion of loci with 0 shared alleles
- `InfType` - Inferred relationship (PO, FS, 2nd, 3rd, UN, Dup/MZ)

**Relationship interpretation:**
| InfType | Kinship | IBS0 | Meaning |
|---------|---------|------|---------|
| PO | ~0.25 | ~0 | Parent-Offspring |
| FS | ~0.25 | >0 | Full Siblings |
| 2nd | ~0.125 | - | Half-sibs, Grandparent |
| UN | <0.05 | - | Unrelated |

---

### PCA Scores CSV

**Usage:** Principal component scores for samples

**Format:** CSV with header
```csv
sample.id,PC1,PC2,PC3,...,PC20
Animal001,-0.0123,0.0456,-0.0789,...,0.0012
Animal002,0.0234,-0.0567,0.0891,...,-0.0023
```

**Columns:**
- `sample.id` - Individual ID
- `PC1` through `PC20` - Principal component scores

---

### Parent-Offspring Summary CSV

**Usage:** Identified parent-offspring pairs

**Format:** CSV with header
```csv
Individual,Parent,InfType,Kinship,IBS0
2190762,1040084,PO,0.2456,0.0012
2210676,2000348,PO,0.2501,0.0008
```

**Columns:**
- `Individual` - Query animal ID
- `Parent` - Identified parent ID
- `InfType` - Relationship type (always "PO")
- `Kinship` - Kinship coefficient (~0.25)
- `IBS0` - IBS0 statistic (~0 for true PO)

---

## Reference File Formats

### SNP Map File

**Usage:** Mapping SNP names to genome coordinates

**Format:** Tab-delimited with header
```
SNP_Name    Chr    Position    Allele1    Allele2
ARS-BFGL-BAC-10172    1    123456    A    G
ARS-BFGL-BAC-10176    1    156789    C    T
```

---

### Alleles Template

**Usage:** Reference allele coding

**Format:** Space-delimited, 3 columns, no header
```
ARS-BFGL-BAC-10172    A    G
ARS-BFGL-BAC-10176    C    T
```

**Purpose:** Ensures consistent allele coding across datasets

---

## Data Conversion Examples

### VCF to PLINK

```bash
plink --vcf input.vcf.gz --make-bed --out output
```

### PLINK to VCF

```bash
plink --bfile input --recode vcf-iid --out output
bcftools view -Oz output.vcf -o output.vcf.gz
tabix -p vcf output.vcf.gz
```

### CSV to PLINK ID List

```bash
# Remove header and extract columns
tail -n +2 input.csv | awk -F',' '{print $1, $2}' > output.list
```

---

## File Naming Conventions

### PLINK Filesrequires prefix:
```
dataset.bed
dataset.bim
dataset.fam
```

### VCF with index:
```
dataset.vcf.gz
dataset.vcf.gz.tbi
```

### Intermediate files:
```
dataset.prune.in
dataset.prune.out
dataset.irem (removed individuals)
dataset.nosex (sex errors)
```

---

## Validation Commands

### Check PLINK Files

```bash
# Verify dimensions
wc -l dataset.bim  # N SNPs
wc -l dataset.fam  # N samples

# Check for errors
plink --bfile dataset --make-just-bim --out check
```

### Check VCF Files

```bash
# Validate format
bcftools view dataset.vcf.gz | head -100

# Count samples
bcftools query -l dataset.vcf.gz | wc -l

# Count variants
bcftools view -H dataset.vcf.gz | wc -l
```

### Check CSV Files

```bash
# View header
head -1 metadata.csv

# Count rows (excluding header)
tail -n +2 metadata.csv | wc -l

# Check for duplicates
awk -F',' 'NR>1{print $2}' metadata.csv | sort | uniq -d
```

---

## Common Format Issues

### Line Endings

**Problem:** Windows (CRLF) vs Unix (LF)

**Detection:**
```bash
file input.txt
# Windows: "ASCII text, with CRLF line terminators"
# Unix: "ASCII text"
```

**Fix:**
```bash
dos2unix input.txt
# OR
sed -i 's/\r$//' input.txt
```

### Character Encoding

**Problem:** Non-ASCII characters

**Detection:**
```bash
file -i input.txt
# Should show: charset=us-ascii or charset=utf-8
```

**Fix:**
```bash
iconv -f ISO-8859-1 -t UTF-8 input.txt > output.txt
```

### Whitespace Issues

**Problem:** Tabs vs spaces, trailing spaces

**Detection:**
```bash
cat -A input.txt  # Shows tabs as ^I, line ends as $
```

**Fix:**
```bash
# Convert tabs to spaces
expand input.txt > output.txt

# Remove trailing whitespace
sed 's/[[:space:]]*$//' input.txt > output.txt
```
