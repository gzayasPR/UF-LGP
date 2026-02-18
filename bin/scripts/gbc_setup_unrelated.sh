#!/bin/bash
#SBATCH --time=48:00:00

#==============================================================================
# ADMIXTURE ANALYSIS - UNRELATED POPULATION SELECTION
#==============================================================================
# Description: Prepares reference and mixed populations for admixture analysis
#              by selecting unrelated individuals and harmonizing SNP sets
#==============================================================================

# Input parameters
proj_env=$1
NAME=$2
IN_GENO=$3
mixedlist=$4
reflist=$5
mind=${6:-0.1}
maf=${7:-0.05}
geno=${8:-0.01}

# Validate required inputs
if [[ -z "$proj_env" || -z "$NAME" || -z "$IN_GENO" || -z "$mixedlist" || -z "$reflist" ]]; then
    echo "ERROR: Missing required arguments"
    echo "Usage: $0 <proj_env> <NAME> <IN_GENO> <mixedlist> <reflist> [mind] [maf] [geno]"
    exit 1
fi

# Source project environment
source "$proj_env"
source "${my_bin}/scripts/shell.functions.sh"
NTHREADS="${SLURM_CPUS_PER_TASK:-12}"

# Load necessary modules
ml plink/1.90b3
ml admixture/1.3.0
ml R/3.6

# Detect input format
PLINK_INPUT=$(detect_genotype_format "$IN_GENO")
echo $PLINK_INPUT
#==============================================================================
# SETUP DIRECTORIES
#==============================================================================
out_dir="${my_results}/${NAME}"
admix_dir="${out_dir}/Mixed.population"
ref_dir="${out_dir}/Reference.Population"
admixture_dir="${out_dir}/ADMIXTURE"
RSCRIPTS="${my_bin}/scripts"

mkdir -p "$admix_dir" "$ref_dir" "$admixture_dir" "${admixture_dir}/Supervised_Unrelated"

echo "=============================================="
echo "Running analysis for: $NAME"
echo "Input genotype: $IN_GENO"
echo "QC parameters: mind=${mind}, maf=${maf}, geno=${geno}"
echo "=============================================="

#==============================================================================
# STEP 1: PREPARE REFERENCE POPULATION
#==============================================================================
echo "[Step 1/5] Processing reference population..."
cd "$ref_dir"

# Create reference sample list
cat "$reflist" > "${ref_dir}/all.list.ID"

# QC and LD pruning for reference population
plink --cow \
    $PLINK_INPUT \
    --keep "${ref_dir}/all.list.ID" \
    --chr 1-29 \
    --mind ${mind} \
    --maf ${maf} \
    --geno ${geno} \
    --indep-pairwise 50 10 0.2 \
    --make-bed \
    --out "${ref_dir}/REF.1" \
    --threads ${NTHREADS}

# Extract pruned SNPs
plink --cow -bfile "${ref_dir}/REF.1" \
    --extract "${ref_dir}/REF.1.prune.in" \
    --make-bed \
    --out "${ref_dir}/REF" \
    --threads ${NTHREADS}

# Create allele reference template (SNP, A1, A2)
awk '{print $2, $5, $6}' "${ref_dir}/REF.bim" > "${ref_dir}/alleles.template"

echo "Reference population: $(wc -l < ${ref_dir}/REF.fam) individuals, $(wc -l < ${ref_dir}/REF.bim) SNPs"

#==============================================================================
# STEP 2: PREPARE MIXED POPULATION ON SAME SNP SET
#==============================================================================
echo "[Step 2/5] Processing mixed population..."
cd "$admix_dir"

# Extract mixed population with same SNP set and allele coding
plink --cow \
     $PLINK_INPUT \
    --keep "$mixedlist" \
    --chr 1-29 \
    --mind ${mind} \
    --geno ${geno} \
    --extract "${ref_dir}/REF.1.prune.in" \
    --a1-allele "${ref_dir}/alleles.template" 3 1 \
    --make-bed \
    --out "${NAME}" \
    --threads ${NTHREADS}

# Write SNP list for mixed population
plink --cow -bfile "${NAME}" \
    --write-snplist \
    --out "${NAME}"

echo "Mixed population: $(wc -l < ${admix_dir}/${NAME}.fam) individuals, $(wc -l < ${admix_dir}/${NAME}.bim) SNPs"

#==============================================================================
# STEP 3: RUN PCA/KING/PC-AiR AND SELECT UNRELATED INDIVIDUALS
#==============================================================================
echo "[Step 3/5] Running relatedness analysis and selecting unrelated individuals..."
cd "$ref_dir"

# Run R script for unrelated individual selection
Rscript "${RSCRIPTS}/Selection_of_Unrelated_script.R" "REF"


#==============================================================================
# STEP 4: BUILD FINAL UNRELATED REFERENCE MATCHED TO MIXED POPULATION
#==============================================================================
echo "[Step 4/5] Building final unrelated reference dataset..."
cd "$ref_dir"

# Create final reference with unrelated individuals, matched SNP order/alleles
plink --cow -bfile "REF" \
    --keep "${ref_dir}/Unrelated_list.txt" \
    --extract "${admix_dir}/${NAME}.snplist" \
    --a1-allele "${ref_dir}/alleles.template" 3 1 \
    --make-bed \
    --out "REF.unrelated" \
    --threads ${NTHREADS}

#==============================================================================
# STEP 5: SANITY CHECKS
#==============================================================================
echo "[Step 5/5] Running sanity checks..."

# Check SNP counts
ref_snps=$(wc -l < "${ref_dir}/REF.unrelated.bim")
mixed_snps=$(wc -l < "${admix_dir}/${NAME}.bim")

echo "REF.unrelated SNPs: ${ref_snps}"
echo "Mixed population SNPs: ${mixed_snps}"

if [[ $ref_snps -ne $mixed_snps ]]; then
    echo "ERROR: SNP count mismatch!"
    exit 1
fi

# Check SNP order and positions
paste "${ref_dir}/REF.unrelated.bim" "${admix_dir}/${NAME}.bim" | \
    awk '{if ($1!=$7 || $4!=$10 || $2!=$8) {
        print "ERROR: SNP mismatch at line", NR;
        print "REF:", $1, $2, $4;
        print "MIXED:", $7, $8, $10;
        exit 1;
    }}'

if [[ $? -ne 0 ]]; then
    echo "ERROR: SNP order/position mismatch detected!"
    exit 1
fi

# Check allele concordance
paste "${ref_dir}/REF.unrelated.bim" "${admix_dir}/${NAME}.bim" | \
    awk '{if ($5!=$11 || $6!=$12) {
        print "ERROR: Allele mismatch at line", NR;
        print "REF:", $2, $5, $6;
        print "MIXED:", $8, $11, $12;
        exit 1;
    }}'

if [[ $? -ne 0 ]]; then
    echo "ERROR: Allele coding mismatch detected!"
    exit 1
fi

#==============================================================================
# COMPLETION
#==============================================================================
echo "=============================================="
echo "✓ All sanity checks passed!"
echo "=============================================="
echo "Final datasets ready for ADMIXTURE analysis:"
echo "  Reference: ${ref_dir}/REF.unrelated"
echo "  Mixed:     ${admix_dir}/${NAME}"
echo "=============================================="

exit 0