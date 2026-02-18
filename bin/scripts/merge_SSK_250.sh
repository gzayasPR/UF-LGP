#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --output=SSK_250K_%j.log
#SBATCH --error=SSK_250K_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8gb
# -------------------- Config --------------------


proj_env=$1
NAME=$2
REF=$3
vcf_ssk=$4
vcf_250k=$5

source $proj_env
results="${my_results}"
out_dir="${results}/${NAME}/merge_SSK_250K"

work="${out_dir}"
mkdir -p "${work}"

ml bcftools/1.22
# Build regions string "1,2,3,...,29" for autosomes
regions="$(seq -s, 1 29)"

echo "INFO: Ensuring BGZF compression for inputs…"
# Re-BGZF to guarantee tabix-compat compression (handles .vcf or non-BGZF .vcf.gz)
ssk_bgz="${work}/SSK.input.bgzf.vcf.gz"
k250_bgz="${work}/250K.input.bgzf.vcf.gz"

bcftools view -Ou "${vcf_ssk}" \
  | bcftools view -Oz --threads "${NTHREADS}" -o "${ssk_bgz}"
tabix -f -p vcf "${ssk_bgz}"

bcftools view -Ou "${vcf_250k}" \
  | bcftools view -Oz --threads "${NTHREADS}" -o "${k250_bgz}"
tabix -f -p vcf "${k250_bgz}"
# -------------------- Find COMMON SNP sites, extract, then merge --------------------

ssk_norm_id="${work}/SSK.norm.id.autosomes.vcf.gz"
k250_norm_id="${work}/250K.norm.id.autosomes.vcf.gz"

common_sites="${work}/SSK_250K.common_snp_sites.txt"

echo "INFO: Finding SNP sites present in BOTH VCFs (intersection of sites)…"
# -Oz: bgzipped VCF output
# -n=2: require presence in both files (2 inputs)
# -v snps: keep SNPs only
bcftools isec \
  -n=2 -c all \
  -Oz -o "${common_sites}" \
  "${ssk_norm_id}" "${k250_norm_id}"


echo "INFO: Extracting common SNP sites from each panel…"
ssk_common="${work}/SSK.common_snps.vcf.gz"
k250_common="${work}/250K.common_snps.vcf.gz"

# Use the intersected VCF as the sites reference. This avoids ID mismatches and handles multiallelic splits cleanly.
bcftools view \
  --threads "${NTHREADS}" \
  -T "${common_sites}" \
  -Oz -o "${ssk_common}" \
  "${ssk_norm_id}"
tabix -f -p vcf "${ssk_common}"

bcftools view \
  --threads "${NTHREADS}" \
  -T "${common_sites}" \
  -Oz -o "${k250_common}" \
  "${k250_norm_id}"
tabix -f -p vcf "${k250_common}"

echo "INFO: Merging the two VCFs on COMMON SNP sites only (samples union; sites shared)…"
merged_union="${out_dir}/SkimSeek_250K.commonSites.merged.vcf.gz"
bcftools merge \
  --threads "${NTHREADS}" \
  --force-samples \
  -Oz -o "${merged_union}" \
  "${ssk_common}" "${k250_common}"

tabix -f -p vcf "${merged_union}"
echo "INFO: Done. Output: ${merged_union}"
