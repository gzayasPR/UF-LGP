# Troubleshooting Guide

Common issues and solutions for UF-LGP pipelines.

## General Issues

### SLURM Job Failures

**Issue:** Job fails immediately or doesn't start

```bash
# Check job status
squeue -u $USER

# View job details
scontrol show job [JOBID]

# Check error log
tail bash_out/[PIPELINE]/[NAME]_[JOBID].log
```

**Common causes:**
- Invalid account name
- Insufficient resources requested
- Module not available
- Syntax error in script

### Module Loading Errors

**Issue:** `module: command not found` or `module not available`

**Solutions:**
```bash
# Check available modules
module avail

# Load specific version
module load plink/1.90b3.39

# Check loaded modules
module list
```

### Path Issues

**Issue:** `No such file or directory`

**Solutions:**
```bash
# Verify environment setup
source bin/project_env.sh
echo $my_bin
echo $my_metadata
echo $my_results

# Check file permissions
ls -la [filepath]

# Use absolute paths if relative paths fail
```

## Pipeline-Specific Issues

### Illumina Processing

#### Issue: Wrong Column Count in FinalReport

**Symptoms:** Pipeline reports column count != 11

**Solution:**
```bash
# Check column count
head -20 [FinalReport.txt] | awk '{print NF}' | sort -u

# If 12 columns, likely extra space in Sample_ID
# Fix with sed (adjust pattern to your case)
sed -i 's/Sample ID/Sample_ID/g' [FinalReport.txt]
```

#### Issue: PLINK Conversion Fails

**Symptoms:** Empty .bed file or PLINK errors

**Solution:**
```bash
# Check LGEN file format
head [NAME].lgen
# Should be: FID IID SNP Allele1 Allele2

# Verify no missing values where not expected
grep -c "0 0" [NAME].lgen  # Count missing genotypes
```

#### Issue: Coordinate Update Fails

**Symptoms:** SNPs remain unmapped after --update-chr

**Solution:**
```bash
# Verify SNP names match between files
head [NAME].bim
head ARS.UF.Chr.Update.txt

# Check for name format differences (spaces, case, etc.)
```

### Skim-Seek VCF Merging

#### Issue: bcftools setGT Plugin Not Found

**Symptoms:** `ERROR: plugin 'setGT' not available`

**Solution:**
```bash
# Option 1: Load bcftools with plugins
ml bcftools/1.22  # or newer with plugins

# Option 2: Skip filtering (keep all genotypes)
# Edit skim_seek.sh: KEEP_LOWCONF=1
```

#### Issue: Memory Exhaustion During Merge

**Symptoms:** Job killed, "out of memory" errors

**Solution:**
```bash
# Increase memory allocation
#SBATCH --mem-per-cpu=16gb  # up from 12gb

# OR reduce parallel workers
CPUS=20  # down from 30

# OR process in chromosome batches
```

#### Issue: Sample Name Mismatches

**Symptoms:** Duplicate or missing samples after merge

**Solution:**
```bash
# Check metadata CSV for duplicates
awk -F',' 'NR>1{print $2}' [meta.csv] | sort | uniq -d

# Verify SkimSeek_ID matches VCF sample name
bcftools query -l [sample.vcf.gz]
```

#### Issue: Contig Name Inconsistencies

**Symptoms:** Missing chromosomes in merged VCF

**Solution:**
```bash
# Check contig names in VCF headers
bcftools view -h [sample.vcf.gz] | grep "^##contig"

# Ensure consistent naming (all "1" or all "chr1", not mixed)
# May need to rename contigs:
bcftools annotate --rename-chrs [chr_name_conv.txt] [input.vcf.gz]
```

### Breed Composition

#### Issue: ADMIXTURE Won't Converge

**Symptoms:** Runs >100 iterations, doesn't finish

**Solutions:**
```bash
# Increase LD pruning stringency
# Edit gbc_setup_unrelated.sh:
--indep-pairwise 50 10 0.1  # stricter (was 0.2)

# Remove high-missingness SNPs
--geno 0.05  # stricter (was 0.1)

# Verify reference populations are distinct
# Check for sample swaps or mislabeling
```

#### Issue: Reference Animals Show Admixture

**Symptoms:** REF.unrelated.2.Q shows Q values far from 0 or 1

**Solutions:**
```bash
# Check which animals are affected
paste REF.unrelated.fam REF.unrelated.2.Q | awk '$3<0.9 && $4<0.9 {print}'

# Verify these are truly purebred
# Remove questionable individuals from reflist
```

#### Issue: PC-AiR Fails to Find Unrelated Set

**Symptoms:** Very few unrelated individuals selected

**Solutions:**
```R
# Adjust kinship threshold in Selection_of_Unrelated_script.R
# In pcair() call:
kin.thresh = 0.05  # More permissive (default 2^(-9/2) ≈ 0.044)
div.thresh = -0.05  # More permissive (default -2^(-9/2))
```

#### Issue: R Package Missing

**Symptoms:** `Error: package 'SNPRelate' not found`

**Solutions:**
```r
# Install missing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SNPRelate", "GENESIS", "gdsfmt"))
install.packages(c("ggplot2", "dplyr", "tidyr"))
```

### Parentage Verification

#### Issue: No PO Relationships Found

**Symptoms:** Empty parent_offspring results

**Solutions:**
```bash
# Verify query IDs in genotype data
awk '{print $2}' king_input.fam | sort > fam_ids.txt
sort check.parentage.ID > query_ids.txt
comm -12 fam_ids.txt query_ids.txt  # Should show matching IDs

# Check reference population size
wc -l [reflist]  # Need potential parents

# Verify KING ran successfully
ls -lh king.kin  # Should be non-empty
```

#### Issue: Unexpected Kinship Values

**Symptoms:** Kinship far from 0.25 for "PO"

**Solutions:**
```bash
# Check genotype quality
plink --bfile king_input --missing --out qc_check
# Review qc_check.imiss and qc_check.lmiss

# Increase QC stringency
--mind 0.05  # stricter
--geno 0.005  # stricter

# Verify sample identities (possible sample swap)
```

#### Issue: High IBS0 for Parent-Offspring

**Symptoms:** IBS0 > 0.01 for InfType="PO"

**Solutions:**
```bash
# May not be true parent-offspring (possibly full-sibs)
# Check kinship.kin file directly
awk '$15=="PO" && $8>0.01' king.kin

# Review these manually - may need pedigree correction
```

## Data Format Issues

### VCF Format Errors

**Issue:** bcftools fails with format errors

**Solutions:**
```bash
# Validate VCF
bcftools view [input.vcf.gz] | head -100

# Re-compress with proper BGZF
bcftools view -Oz [input.vcf] -o [output.vcf.gz]
tabix -p vcf [output.vcf.gz]

# Check for malformed lines
zcat [input.vcf.gz] | grep -v "^#" | awk 'NF!=9+nsamples'
```

### PLINK Format Mismatches

**Issue:** .bed/.bim/.fam files inconsistent

**Solutions:**
```bash
# Verify dimensions match
wc -l [prefix].bim  # N SNPs
wc -l [prefix].fam  # N samples

# .bed file size should be roughly: N_SNPs * N_Samples / 4 bytes
# Check:
ls -lh [prefix].bed
```

### ID Format Problems

**Issue:** FID/IID not matching expected format

**Solutions:**
```bash
# PLINK expects: FID IID (space or tab separated)
# Check format
head [id_list.txt]

# Fix if comma-separated
sed 's/,/\t/g' [id_list.txt] > [id_list_fixed.txt]

# Remove quotes if present
sed 's/"//g' [id_list.txt] > [id_list_fixed.txt]
```

## Performance Issues

### Slow Runtime

**Issue:** Job takes much longer than expected

**Solutions:**
```bash
# Check if using parallel processing
echo $SLURM_CPUS_PER_TASK

# Monitor CPU usage
sstat -j [JOBID] --format=AveCPU,MaxRSS

# Optimize LD pruning (fewer SNPs = faster)
--indep-pairwise 50 10 0.2

# Use faster algorithms where available
--threads [N]  # for tools that support it
```

### Disk Space Issues

**Issue:** No space left on device

**Solutions:**
```bash
# Check disk usage
df -h $my_results
du -sh $my_results/*

# Clean up intermediate files
rm -rf results/[NAME]/Make_PLINK/
rm -rf results/[NAME]/*/temp*

# Use scratch space for large temp files
scratch="/scratch/local/${SLURM_JOB_ID}"
```

### Memory Issues

**Issue:** Job killed due to memory

**Solutions:**
```bash
# Check memory usage in logs
grep -i "memory" bash_out/[PIPELINE]/[NAME]_*.log

# Increase allocation
#SBATCH --mem-per-cpu=16gb

# Process in smaller batches
# Or reduce parallel workers
```

## Getting Help

If you can't resolve an issue:

1. **Check logs:**
   ```bash
   tail -100 bash_out/[PIPELINE]/[NAME]_[JOBID].log
   ```

2. **Verify inputs:**
   ```bash
   head [input_file]
   wc -l [input_file]
   ```

3. **Test with small dataset:**
   - Subset to chromosome 29 (smallest)
   - Use 10-100 samples
   - Verify pipeline works before scaling up

4. **Collect information:**
   - Error messages (exact text)
   - Command that failed
   - Input file formats
   - Software versions: `plink --version`, `bcftools --version`

5. **Create issue on GitHub:**
   - Include error log
   - Describe what you've tried
   - Provide example data if possible

## Preventive Measures

### Before Running Pipelines

1. **Validate inputs:**
   ```bash
   # Check files exist
   ls -lh [input_files]
   
   # Check formats
   head [input_files]
   ```

2. **Test with subset:**
   ```bash
   # Extract chr29 only
   plink --bfile [input] --chr 29 --make-bed --out test_chr29
   ```

3. **Verify environment:**
   ```bash
   # Check modules
   ml plink bcftools R
   
   # Verify paths
   echo $my_bin $my_metadata $my_results
   ```

4. **Review configuration:**
   ```bash
   # Check SLURM parameters are appropriate
   # Verify file paths are correct
   ```

### During Pipeline Execution

1. **Monitor jobs:**
   ```bash
   watch squeue -u $USER
   ```

2. **Check logs periodically:**
   ```bash
   tail -f bash_out/[PIPELINE]/[NAME]_*.log
   ```

3. **Verify intermediate outputs:**
   ```bash
   ls -lh results/[NAME]/
   ```

### After Pipeline Completion

1. **Validate outputs:**
   ```bash
   # Check dimensions
   wc -l results/[NAME]/[output_files]
   
   # Check quality
   head results/[NAME]/[output_files]
   ```

2. **Archive if successful:**
   ```bash
   # Save important outputs
   # Delete intermediate files
   ```

3. **Document parameters used:**
   ```bash
   # Keep record of settings for reproducibility
   ```
