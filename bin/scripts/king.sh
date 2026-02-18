#!/bin/bash
#SBATCH --time=48:00:00

# Input parameters
proj_env=$1
NAME=$2
IN_GENO=$3
reflist=$4
INDVI_LIST=$5
mind=${6:-0.1}
maf=${7:-0.05}
geno=${8:-0.01}

# Validate required inputs
if [[ -z "$proj_env" || -z "$NAME" || -z "$IN_GENO" || -z "$reflist" || -z "$INDVI_LIST" ]]; then
    echo "ERROR: Missing required arguments"
    echo "Usage: $0 <proj_env> <NAME> <IN_GENO> <reflist> <INDVI_LIST> [mind] [maf] [geno]"
    exit 1
fi

# Validate individual list file exists
if [[ ! -f "$INDVI_LIST" ]]; then
    echo "ERROR: Individual list file not found: $INDVI_LIST"
    exit 1
fi

# Source project environment
source "$proj_env"
NTHREADS="${SLURM_CPUS_PER_TASK:-12}"

# Setup working directory
WDIR="${my_results}/${NAME}"
mkdir -p "$WDIR"
cd "$WDIR"

echo "Starting parentage assignment for individuals in: ${INDVI_LIST}"
echo "Working directory: ${WDIR}"

# Load PLINK module
ml plink/1.90b3.39

# Step 1: Process and prune the data
echo "Processing and pruning genotype data..."
plink --cow --bfile "${IN_GENO}" --keep "${reflist}" \
      --chr 1-29 --maf ${maf} --geno ${geno} --mind ${mind} --hwe 0.00001 \
      --indep-pairwise 5000 10 0.5 --make-bed --out Temp

# Step 2: Extract pruned SNPs
echo "Extracting pruned SNPs..."
plink --cow --bfile Temp --extract Temp.prune.in --make-bed --out king_input

# Step 3: Run KING for relatedness
echo "Running KING relatedness analysis..."
ml king/2.3.0
king -b king_input.bed --related --sexchr 30

# Step 4: Extract parentage results for each individual in the list
echo "Extracting parentage assignments for individuals..."
mkdir -p parentage_results

# Get header from king.kin
header=$(head -n1 king.kin)

# Process each individual in the list
while IFS= read -r INDVI; do
    # Skip empty lines
    [[ -z "$INDVI" ]] && continue
    
    echo "Processing individual: ${INDVI}"
    
    # Extract results for this individual
    grep "${INDVI}" king.kin > "parentage_results/${INDVI}.kin"
    
    # Add header if results were found
    if [[ -s "parentage_results/${INDVI}.kin" ]]; then
        sed -i "1s/^/${header}\n/" "parentage_results/${INDVI}.kin"
        echo "  Found relatedness data for ${INDVI}"
    else
        echo "  WARNING: No relatedness data found for ${INDVI}"
        rm -f "parentage_results/${INDVI}.kin"
    fi
done < "$INDVI_LIST"

# Create a combined results file with all individuals
echo "Creating combined results file..."
echo "$header" > parentage_results/all_individuals.kin
while IFS= read -r INDVI; do
    [[ -z "$INDVI" ]] && continue
    grep "${INDVI}" king.kin >> parentage_results/all_individuals.kin
done < "$INDVI_LIST"

# Step 5: Extract parent-offspring relationships
echo ""
echo "=== IDENTIFYING PARENT-OFFSPRING RELATIONSHIPS ==="
echo ""

# Create parent-offspring results directory
mkdir -p parentage_results/parent_offspring

# Header for parent-offspring summary
echo "Individual,Parent/Offspring,InfType,Kinship,IBS0" > parent_offspring_summary.csv
# Process each individual to find parent-offspring relationships
while IFS= read -r INDVI; do
    [[ -z "$INDVI" ]] && continue
    
    if [[ -f "parentage_results/${INDVI}.kin" ]]; then
        echo "Checking for parents of: ${INDVI}"
        
        # Extract parent-offspring relationships (InfType = PO)
        # KING uses "PO" to indicate parent-offspring relationships
        # InfType is column 15
        awk 'NR>1 && $15=="PO"' "parentage_results/${INDVI}.kin" > "parentage_results/parent_offspring/${INDVI}_parents.txt"
        
        if [[ -s "parentage_results/parent_offspring/${INDVI}_parents.txt" ]]; then
            # Add header to parent file
            echo "$header" | cat - "parentage_results/parent_offspring/${INDVI}_parents.txt" > temp_file
            mv temp_file "parentage_results/parent_offspring/${INDVI}_parents.txt"
            
            # Count parents found
            parent_count=$(wc -l < "parentage_results/parent_offspring/${INDVI}_parents.txt")
            parent_count=$((parent_count - 1))  # Subtract header
            echo "  ✓ Found ${parent_count} parent(s) for ${INDVI}"
            
            # Extract parent IDs and add to summary
            # Columns: $2=ID1, $3=ID2, $8=IBS0, $11=Kinship, $15=InfType
            awk -v indv="$INDVI" 'NR>1 {
                # Determine which ID is the parent (the one that is NOT the query individual)
                if ($2 == indv) {
                    parent_id = $3
                } else {
                    parent_id = $2
                }
                printf "%s,%s,%s,%s,%s\n", indv, parent_id, $15, $11, $8
            }' "parentage_results/parent_offspring/${INDVI}_parents.txt" >> parent_offspring_summary.csv
            
            # Display to console
            awk -v indv="$INDVI" 'NR>1 {
                if ($2 == indv) parent_id = $3; else parent_id = $2
                printf "    Parent: %s (Kinship: %s, IBS0: %s)\n", parent_id, $11, $8
            }' "parentage_results/parent_offspring/${INDVI}_parents.txt"
        else
            echo "  ✗ No parent-offspring relationships found for ${INDVI}"
            rm -f "parentage_results/parent_offspring/${INDVI}_parents.txt"
        fi
    fi
done < "$INDVI_LIST"

# Create a summary report
echo ""
echo "=== PARENTAGE ASSIGNMENT SUMMARY ==="
cat parentage_results/parent_offspring_summary.csv

# Step 6: Cleanup intermediate files
echo ""
echo "Cleaning up intermediate files..."
rm -f Temp.bed Temp.bim Temp.fam Temp.log
rm -f Temp.prune.in Temp.prune.out
rm -f king_input.bed king_input.bim king_input.fam
rm -f king.kin king.kin0
rm -f *.nosex *.log
rm -f Temp.irem kingallsegs.txt
echo ""
echo "========================================="
echo "Parentage assignment complete!"
echo "========================================="
echo "Individual results: parentage_results/"
echo "Combined results: parentage_results/all_individuals.kin"
echo "Parent-offspring relationships: parentage_results/parent_offspring/"
echo "Summary CSV: parentage_results/parent_offspring_summary.csv"
echo "========================================="