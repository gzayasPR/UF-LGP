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
# ============================================================
#  Neogen SKIM_seek per-sample VCF merge (confidence-aware)
#
#  Key change vs your original:
#    - Keep both PASS and LOWCONF records during cleaning
#    - Toggle behavior for low-confidence genotypes:
#         KEEP_LOWCONF=1  -> keep LOWCONF genotypes as-called
#         KEEP_LOWCONF=0  -> set low-confidence genotypes to missing (./.)
#            Low-confidence defined as:
#              FILTER="LOWCONF" OR MAX(FMT/GP) < GP_MIN
#
#  Then strip everything down to GT only (plus headers), so merges are stable.
#
#  Metadata CSV columns (header required):
#    path,SkimSeek_ID,ProjectID,IID,NewID
# ============================================================
# -------------------- Paths / Config --------------------
proj_env="/blue/mateescu/gzayas97/UF-LGP/bin/project_env.sh"
NAME="Skim_Seek_2026"
meta_csv="/blue/mateescu/gzayas97/For_UF/Skim_Seq/data/UF_SkimSeek/SSK_IDs_2026.csv"

account="mateescu"
CPUS=30
mem_CPU=12Gb
max_time=120:00:00

source $proj_env

# Step 1) Merge VCF Files
# -------------------- Confidence toggles --------------------
# 1 = keep LOWCONF and low-GP genotypes as called
# 0 = set them to missing GT=./.
KEEP_LOWCONF=0
# Used only when KEEP_LOWCONF=0
GP_MIN=0.95
SSK_log=${my_bash}/Skim_Seek/
mkdir -p $SSK_log
sbatch --account $account  \
      --job-name="SSK_${NAME}" \
      --cpus-per-task $CPUS \
      --mem-per-cpu $mem_CPU  \
      --time=${max_time} \
      --output="${SSK_log}/${NAME}%j.log"  \
      --error="${SSK_log}/${NAME}%j.log"  \
      ${my_bin}/scripts/merge_skim_seek.sh \
      $proj_env \
      $NAME \
      $meta_csv \
      $KEEP_LOWCONF \
      $GP_MIN


