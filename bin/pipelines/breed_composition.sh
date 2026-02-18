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
NAME="GBC_Febuary_2026"
IN_GENO="/blue/mateescu/gzayas97/For_UF/Skim_Seq/data/250K/UF_250K"
mixedlist=/blue/mateescu/gzayas97/UF-LGP/metadata/IDs/admix.ID
reflist=/blue/mateescu/gzayas97/UF-LGP/metadata/IDs/MAB.ID
account="mateescu"
CPUS=3
mem_CPU=12Gb
max_time=120:00:00
maf=0.01
mind=0.1
geno=0.1


source $proj_env

GBC_log=${my_bash}/${NAME}/
mkdir -p $GBC_log

sbatch --account $account  \
      --job-name="Est_Unrel_${NAME}" \
      --cpus-per-task $CPUS \
      --mem-per-cpu $mem_CPU  \
      --time=${max_time} \
      --output="${GBC_log}/Est_Unrel_${NAME}%j.log"  \
      --error="${GBC_log}/Est_Unrel_${NAME}%j.log"  \
      ${my_bin}/scripts/gbc_setup_unrelated.sh \
      $proj_env \
      $NAME \
      $IN_GENO \
      $mixedlist \
      $reflist \
      $mind \
      $maf \
      $geno

sbatch --account $account  \
      --job-name="ADMIX_${NAME}" \
      --cpus-per-task $CPUS \
      --mem-per-cpu $mem_CPU  \
      --time=${max_time} \
      --output="${GBC_log}/ADMIX_${NAME}%j.log"  \
      --error="${GBC_log}/ADMIX_${NAME}%j.log"  \
      ${my_bin}/scripts/gbc_admixture.sh \
      $proj_env \
      $NAME \
      $IN_GENO \
      $mixedlist \
      $reflist \
      $mind \
      $maf \
      $geno
