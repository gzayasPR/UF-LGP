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
REF="/blue/mateescu/gzayas97/genome_refs/bos_taurus/ensembl/ARS-UCD1.2/109/cdna_no_contigs/Bos_taurus.ARS-UCD1.2.cdna.all.fa.no_contigs.fa"
vcf_ssk="/blue/mateescu/gzayas97/UF-LGP/results/Skim_Seek_2026/Skim_Seek_2026_merged_output.vcf.gz"
vcf_250k=

account="mateescu"
CPUS=5
mem_CPU=12Gb
max_time=120:00:00

source $proj_env

SSK_log=${my_bash}/Merge_250K_SSK/
mkdir -p $SSK_log
sbatch --account $account  \
      --job-name="SSK_${NAME}" \
      --cpus-per-task $CPUS \
      --mem-per-cpu $mem_CPU  \
      --time=${max_time} \
      --output="${SSK_log}/${NAME}%j.log"  \
      --error="${SSK_log}/${NAME}%j.log"  \
      ${my_bin}/scripts/merge_SSK_250.sh \
      $proj_env \
      $NAME \
      $REF \
      $vcf_ssk \
      $vcf_250k

