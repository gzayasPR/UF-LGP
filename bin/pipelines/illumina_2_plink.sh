#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Illumina.Prep
#SBATCH --output=Illumina.Prep_%j.out
#SBATCH --error=errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

proj_env="/blue/mateescu/gzayas97/UF-LGP/bin/project_env.sh"
NAME="250K_2026"
Illumina_meta_csv="/blue/mateescu/gzayas97/UF-LGP/metadata/250K_finalreports.csv"
SNP_MAP="/blue/mateescu/gzayas97/UF-LGP/metadata/Illumina_to_Plink/SNP_Map.txt"
UPDATE_IDS="/blue/mateescu/gzayas97/UF-LGP/metadata/Illumina_to_Plink/Feb_2024.txt"
CHR_UPD="/blue/mateescu/gzayas97/UF-LGP/metadata/Illumina_to_Plink/ARS.UF.Chr.Update.txt"
BP_UPD="/blue/mateescu/gzayas97/UF-LGP/metadata/Illumina_to_Plink/ARS.BP.2024.txt"
INDELS="/blue/mateescu/gzayas97/UF-LGP/metadata/Illumina_to_Plink/INDELS.ID"
ALLELE_POS=4 # 4 is the ATCG and 8 is the Illumina A and B allele format
REF_GENOME="/blue/mateescu/gzayas97/genome_refs/bos_taurus/ensembl/ARS-UCD1.2/109/cdna_no_contigs/Bos_taurus.ARS-UCD1.2.cdna.all.fa.no_contigs.fa"

account="mateescu"
CPUS=10
mem_CPU=12Gb
max_time=120:00:00

source $proj_env
Ill_2_PLK_log=${my_bash}/Illumina_2_PLINK/
mkdir -p $Ill_2_PLK_log

sbatch --account $account  \
      --job-name="Il_prep_${NAME}" \
      --cpus-per-task $CPUS \
      --mem-per-cpu $mem_CPU  \
      --time=${max_time} \
      --output="${Ill_2_PLK_log}/PrepIll_${NAME}%j.log"  \
      --error="${Ill_2_PLK_log}/PrepIll_${NAME}%j.log"  \
      ${my_bin}/scripts/illumina_prep.sh \
      $proj_env \
      $NAME \
      $Illumina_meta_csv

sbatch --account $account  \
      --job-name="I2L_${NAME}" \
      --cpus-per-task $CPUS \
      --mem-per-cpu $mem_CPU  \
      --time=${max_time} \
      --output="${Ill_2_PLK_log}/makePLINK_${NAME}%j.log"  \
      --error="${Ill_2_PLK_log}/makePLINK_${NAME}%j.log"  \
      ${my_bin}/scripts/make_ill_plink.sh \
      $proj_env \
      $NAME \
      $SNP_MAP \
      $UPDATE_IDS \
      $CHR_UPD \
      $BP_UPD \
      $INDELS \
      $ALLELE_POS \
      $REF_GENOME