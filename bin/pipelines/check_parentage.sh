#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --output=Making_VCF_%j.log
#SBATCH --error=Making_VCF_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb

######################## UF Livestock Genomics Pipelines (UF-LGP) ##############3333
################################################################################
# PARENTAGE ASSIGNMENT PIPELINE USING KING
################################################################################
# 
# DESCRIPTION:
#   This script performs parentage assignment analysis using KING relatedness
#   inference. It processes a list of individuals, filters genotype data,
#   runs KING relatedness analysis, and identifies parent-offspring pairs.
#
# USAGE:
#   sbatch script.sh <proj_env> <NAME> <IN_GENO> <reflist> <INDVI_LIST> [mind] [maf] [geno]
#
# REQUIRED ARGUMENTS:
#   proj_env    : Path to project environment file to source
#   NAME        : Output directory name
#   IN_GENO     : Path to input PLINK binary fileset (without extension)
#   reflist     : File with reference individuals to keep (FID IID format)
#   INDVI_LIST  : File with individual IDs to analyze (one per line)
#
# OPTIONAL ARGUMENTS:
#   mind        : Maximum per-individual missing rate (default: 0.1)
#   maf         : Minor allele frequency threshold (default: 0.05)
#   geno        : Maximum per-SNP missing rate (default: 0.01)
#
proj_env="/blue/mateescu/gzayas97/UF-LGP/bin/project_env.sh"
NAME="Parentage_2026"
IN_GENO="/blue/mateescu/gzayas97/For_UF/Skim_Seq/data/250K/UF_250K"
reflist=/blue/mateescu/gzayas97/UF-LGP/metadata/IDs/MAB.ID
INDVI_LIST=/blue/mateescu/gzayas97/UF-LGP/metadata/check.parentage.ID
account="mateescu"
CPUS=3
mem_CPU=12Gb
max_time=120:00:00
maf=0.01
mind=0.1
geno=0.1


source $proj_env

PARENT_log=${my_bash}/${NAME}/
mkdir -p $PARENT_log

sbatch --account $account  \
      --job-name="KING_${NAME}" \
      --cpus-per-task $CPUS \
      --mem-per-cpu $mem_CPU  \
      --time=${max_time} \
      --output="${PARENT_log}/KING_${NAME}%j.log"  \
      --error="${PARENT_log}/KING_${NAME}%j.log"  \
      ${my_bin}/scripts/king.sh \
      $proj_env \
      $NAME \
      $IN_GENO \
      $reflist \
      $INDVI_LIST \
      $mind \
      $maf \
      $geno

################################################################################
# OUTPUT STRUCTURE:
################################################################################
#
# ${my_results}/${NAME}/
# ├── parentage_results/                    # Main results folder
# │   ├── 2190762.kin                      # Individual KING output files
# │   ├── 2210676.kin                      # (all relatedness for each animal)
# │   ├── 2210764.kin
# │   ├── all_individuals.kin              # Combined KING output for all animals
# │   └── parent_offspring/                # Parent-offspring specific results
# │       ├── 2190762_parents.txt          # PO relationships for each animal
# │       ├── 2210676_parents.txt
# │       └── 2210764_parents.txt
# ├── parent_offspring_summary.csv         # Summary CSV with all PO pairs
# └── Temp.irem                            # Animals excluded due to missingness
#
################################################################################
# KING .kin FILE COLUMN DESCRIPTION:
################################################################################
#
# Column | Field      | Description
# -------|------------|----------------------------------------------------------
#   1    | FID        | Family ID - identifier for the family/population group
#   2    | ID1        | Individual ID for the first person in the pair
#   3    | ID2        | Individual ID for the second person in the pair
#   4    | N_SNP      | Number of SNPs used in the relationship analysis
#   5    | Z0         | Probability that the pair shares 0 alleles IBD
#   6    | Phi        | Kinship coefficient (alternative calculation)
#   7    | HetHet     | Proportion of SNPs where both individuals are heterozygous
#   8    | IBS0       | Proportion of SNPs with zero alleles shared IBS
#   9    | HetConc    | Heterozygous concordance rate
#  10    | HomIBS0    | Homozygous IBS0 - used for relationship inference
#  11    | Kinship    | Kinship coefficient (primary relatedness measure)
#  12    | IBD1Seg    | Number of IBD1 segments (sharing 1 allele)
#  13    | IBD2Seg    | Number of IBD2 segments (sharing 2 alleles)
#  14    | PropIBD    | Proportion of genome shared IBD
#  15    | InfType    | Inferred relationship type (see below)
#  16    | Error      | Error flag (0 = no error)
#
# InfType Values:
#   PO      = Parent-Offspring
#   FS      = Full-Sibling
#   2nd     = Second-degree relative
#   3rd     = Third-degree relative
#   UN      = Unrelated
#   Dup/MZ  = Duplicate/Monozygotic twins
#
# KEY COLUMNS FOR PARENTAGE ASSIGNMENT:
#   - InfType (column 15): "PO" indicates parent-offspring relationship
#   - Kinship (column 11): ~0.25 for parent-offspring and full siblings
#                          ~0.125 for 2nd degree relatives
#   - IBS0 (column 8):     Should be very low (~0) for parent-offspring pairs
#                          (parents share at least one allele at every locus)
#
# parent_offspring_summary.csv FORMAT:
#   Individual,Parent,InfType,Kinship,IBS0
#   2190762,1040084,PO,0.2456,0.0012
#   2210676,2000348,PO,0.2501,0.0008
#
################################################################################
