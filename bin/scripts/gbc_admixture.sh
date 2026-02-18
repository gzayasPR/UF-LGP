#!/bin/bash
#SBATCH --time=48:00:00

#==============================================================================
# SUPERVISED ADMIXTURE ANALYSIS WITH VISUALIZATION
#==============================================================================
# Description: Runs supervised ADMIXTURE analysis on unrelated populations
#              and generates Q-plots and PCA visualizations
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


proj_env=/blue/mateescu/gzayas97/UF-LGP/bin/project_env.sh 
NAME=GBC_Febuary_2026 
IN_GENO=/blue/mateescu/gzayas97/For_UF/Skim_Seq/data/250K/UF_250K
mixedlist=/blue/mateescu/gzayas97/UF-LGP/metadata/IDs/admix.ID
reflist=/blue/mateescu/gzayas97/UF-LGP/metadata/IDs/MAB.ID
maf=0.01
mind=0.1
geno=0.1


# Validate required inputs
if [[ -z "$proj_env" || -z "$NAME" || -z "$IN_GENO" || -z "$mixedlist" || -z "$reflist" ]]; then
    echo "ERROR: Missing required arguments"
    echo "Usage: $0 <proj_env> <NAME> <IN_GENO> <mixedlist> <reflist> [mind] [maf] [geno]"
    exit 1
fi

# Source project environment
source "$proj_env"
NTHREADS="${SLURM_CPUS_PER_TASK:-12}"

# Load necessary modules
ml plink/1.90b3
ml admixture/1.3.0
ml R/3.6

# Source project environment
source "$proj_env"
source "${my_bin}/scripts/shell.functions.sh"

# Detect input format
PLINK_INPUT=$(detect_genotype_format "$IN_GENO")

#==============================================================================
# SETUP DIRECTORIES
#==============================================================================
out_dir="${my_results}/${NAME}"
admix_dir="${out_dir}/Mixed.population"
ref_dir="${out_dir}/Reference.Population"
admixture_dir="${out_dir}/ADMIXTURE"
supervised_dir="${admixture_dir}/Supervised_Unrelated"
figures_dir="${out_dir}/Figures"
RSCRIPTS="${my_bin}/scripts"

mkdir -p "$admix_dir" "$ref_dir" "$admixture_dir" "$supervised_dir" "$figures_dir"

echo "=============================================="
echo "Running ADMIXTURE analysis for: $NAME"
echo "Input genotype: $IN_GENO"
echo "K = 2 (supervised mode)"
echo "=============================================="

#==============================================================================
# VALIDATE INPUT FILES
#==============================================================================
echo "[Validation] Checking required input files..."

# Check if reference unrelated dataset exists
if [[ ! -f "${ref_dir}/REF.unrelated.bed" ]]; then
    echo "ERROR: Reference unrelated dataset not found: ${ref_dir}/REF.unrelated.bed"
    echo "Please run the population preparation script first."
    exit 1
fi

# Check if mixed population dataset exists
if [[ ! -f "${admix_dir}/${NAME}.bed" ]]; then
    echo "ERROR: Mixed population dataset not found: ${admix_dir}/${NAME}.bed"
    echo "Please run the population preparation script first."
    exit 1
fi

# Check if unrelated list exists
if [[ ! -f "${ref_dir}/Unrelated_list.txt" ]]; then
    echo "WARNING: Unrelated list not found at ${ref_dir}/Unrelated_list.txt"
    echo "Checking alternative location: ${ref_dir}/Unrelated_list.txt"
    if [[ ! -f "${ref_dir}/Unrelated_list.txt" ]]; then
        echo "ERROR: Unrelated list not found in either location."
        exit 1
    fi
    UNRELATED_LIST="${ref_dir}/Unrelated_list.txt"
else
    UNRELATED_LIST="${ref_dir}/Unrelated_list.txt"
fi

# Check if PCA results exist
if [[ ! -f "${ref_dir}/PC_scores_PC1-20.csv" ]]; then
    echo "WARNING: PCA results not found at ${ref_dir}/PC_scores_PC1-20.csv"
    echo "Checking alternative location: ${ref_dir}/PC_scores_PC1-20.csv"
    if [[ ! -f "${ref_dir}/PC_scores_PC1-20.csv" ]]; then
        echo "WARNING: PCA results not found. PCA visualization will be skipped."
        SKIP_PCA=true
    else
        PCA_FILE="${ref_dir}/PC_scores_PC1-20.csv"
        SKIP_PCA=false
    fi
else
    PCA_FILE="${ref_dir}/PC_scores_PC1-20.csv"
    SKIP_PCA=false
fi

echo "✓ Input validation complete"

#==============================================================================
# STEP 1: RUN SUPERVISED ADMIXTURE ON REFERENCE POPULATION
#==============================================================================
echo "[Step 1/4] Running ADMIXTURE on reference (unrelated) population..."
cd "$supervised_dir"

# Run ADMIXTURE on reference to get allele frequency estimates
admixture "${ref_dir}/REF.unrelated.bed" 2 -j${NTHREADS}

if [[ $? -ne 0 ]]; then
    echo "ERROR: ADMIXTURE failed on reference population"
    exit 1
fi

echo "✓ Reference ADMIXTURE complete"

#==============================================================================
# STEP 2: PREPARE SUPERVISED ANALYSIS FOR MIXED POPULATION
#==============================================================================
echo "[Step 2/4] Preparing supervised analysis for mixed population..."

# Copy reference allele frequencies for supervised mode
cp REF.unrelated.2.P ${NAME}.2.P.in

if [[ ! -f "${NAME}.2.P.in" ]]; then
    echo "ERROR: Failed to create supervised input file: ${NAME}.2.P.in"
    exit 1
fi

echo "✓ Supervised mode prepared (using reference allele frequencies)"

#==============================================================================
# STEP 3: RUN SUPERVISED ADMIXTURE ON MIXED POPULATION
#==============================================================================
echo "[Step 3/4] Running supervised ADMIXTURE on mixed population..."

# Run ADMIXTURE in supervised mode (-P flag)
admixture -P "${admix_dir}/${NAME}.bed" 2 -j${NTHREADS}

if [[ $? -ne 0 ]]; then
    echo "ERROR: Supervised ADMIXTURE failed on mixed population"
    exit 1
fi

echo "✓ Supervised ADMIXTURE complete"

#==============================================================================
# STEP 4: COMBINE RESULTS AND GENERATE VISUALIZATIONS
#==============================================================================
echo "[Step 4/4] Combining results and generating visualizations..."

# Combine Q matrix with FAM file for sample identification
Rscript "${RSCRIPTS}/combine_Q_FAM.R" \
    "${admix_dir}/${NAME}.fam" \
    "${supervised_dir}/${NAME}.2.Q" \
    "${supervised_dir}/${NAME}.combined.ancestry.csv"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to combine Q matrix with FAM file"
    exit 1
fi

echo "✓ Ancestry results combined"

# Generate admixture Q-plot
echo "Generating admixture Q-plot..."
Rscript "${RSCRIPTS}/Q.PLOT.R" \
    "${supervised_dir}/${NAME}.combined.ancestry.csv" \
    "${figures_dir}/admixture.png" \
    "${UNRELATED_LIST}"

if [[ $? -ne 0 ]]; then
    echo "WARNING: Q-plot generation failed"
else
    echo "✓ Q-plot saved to: ${figures_dir}/admixture.png"
fi

# Generate PCA plot with admixture overlay
if [[ "$SKIP_PCA" == false ]]; then
    echo "Generating PCA plot with admixture proportions..."
    Rscript "${RSCRIPTS}/PCA.plots.R" \
        "${PCA_FILE}" \
        "${supervised_dir}/${NAME}.combined.ancestry.csv" \
        "${UNRELATED_LIST}" \
        "${figures_dir}/PCA_Admixture.png"
    
    if [[ $? -ne 0 ]]; then
        echo "WARNING: PCA plot generation failed"
    else
        echo "✓ PCA plot saved to: ${figures_dir}/PCA_Admixture.png"
    fi
else
    echo "⊗ PCA visualization skipped (PCA results not found)"
fi

#==============================================================================
# SUMMARY OF RESULTS
#==============================================================================
echo "=============================================="
echo "✓ ADMIXTURE ANALYSIS COMPLETE"
echo "=============================================="
echo "Output files:"
echo "  Ancestry proportions: ${supervised_dir}/${NAME}.combined.ancestry.csv"
echo "  Q matrix:             ${supervised_dir}/${NAME}.2.Q"
echo "  P matrix:             ${supervised_dir}/REF.unrelated.2.P"
echo ""
echo "Visualizations:"
if [[ -f "${figures_dir}/admixture.png" ]]; then
    echo "  ✓ Q-plot:             ${figures_dir}/admixture.png"
else
    echo "  ✗ Q-plot:             FAILED"
fi
if [[ "$SKIP_PCA" == false && -f "${figures_dir}/PCA_Admixture.png" ]]; then
    echo "  ✓ PCA plot:           ${figures_dir}/PCA_Admixture.png"
elif [[ "$SKIP_PCA" == true ]]; then
    echo "  ⊗ PCA plot:           SKIPPED (no PCA data)"
else
    echo "  ✗ PCA plot:           FAILED"
fi
echo "=============================================="

exit 0