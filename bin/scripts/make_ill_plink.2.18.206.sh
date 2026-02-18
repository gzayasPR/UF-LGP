#!/bin/bash
#SBATCH --time=48:00:00

proj_env=$1
NAME=$2
SNP_MAP=$3
UPDATE_IDS=$4
CHR_UPD=$5
BP_UPD=$6
INDELS=$7
ALLELE_POS=$8

source "$proj_env"

out_dir="${my_results}/${NAME}"
prep_dir="${out_dir}/Illumina"
plink_dir="${out_dir}/Make_PLINK"
mkdir -p "$out_dir" "$prep_dir" "$plink_dir"

if [[ "$ALLELE_POS" != "4" && "$ALLELE_POS" != "8" ]]; then
  echo "WARNING: ALLELE_POS is $ALLELE_POS. Expected 4 (real) or 8 (Illumina). Continuing anyway." >&2
fi

[[ -d "$prep_dir" ]] || { echo "ERROR: missing Illumina dir: $prep_dir" >&2; exit 1; }

mkdir -p "$plink_dir"
cd "$plink_dir"

rm -f all.txt Geno.all.txt ID.txt "${NAME}.lgen" Geno.all.txt2 "${NAME}.fam" "${NAME}.map"

date
echo "Start all.txt"
# combine all FinalReports (remove header rows 1-10 in each file)
illumina_files=( "${prep_dir}"/*.txt )
[[ ${#illumina_files[@]} -gt 0 ]] || { echo "ERROR: no Illumina files matched in $prep_dir" >&2; exit 1; }

tail -n +11 -q "${illumina_files[@]}" >> all.txt

echo
echo "Line counts (raw files):"
wc -l "${illumina_files[@]}" || true
echo "Line count (all.txt):"
wc -l all.txt || true
echo

date
echo "Start Geno.all.txt"
awk -F $'\t' -v a="$ALLELE_POS" 'BEGIN{OFS="\t"} {print $1,$2,$3,$a}' all.txt > Geno.all.txt
# Build ID.txt: Date + SampleID (col2)
# safer loop than: for file in `ls ...`
rm -f ID.txt
for file in "${illumina_files[@]}"; do
  b=$(basename "$file")
  DATE=$(echo "$b" | awk -F"BOVF250V1" '/BOVF250V1/{print $2}' | sed 's/[^0-9]*//g')
  awk -v Date="$DATE" -F $'\t' 'BEGIN { OFS="\t" } {print Date,$2}' "$file" | tail -n +11 >> ID.txt
done

echo "Line count (ID.txt):"
wc -l ID.txt || true
echo "Line count (Geno.all.txt):"
wc -l Geno.all.txt || true
echo

date
echo "Start LGEN file"
paste ID.txt <(awk -F $'\t' 'BEGIN { OFS="\t" } {print $1,$3,$4}' Geno.all.txt ) > Geno.all.txt2

head Geno.all.txt2 || true
sed $'s/-\t-/0\t0/g' Geno.all.txt2 > "${NAME}.lgen"

# MAP
cp -f "$SNP_MAP" .
awk -F $'\t' 'BEGIN { OFS="\t" } NR>1 {print $3,$2,0,$4}' "$(basename "$SNP_MAP")" > "${NAME}.map"

# FAM
date
echo "Build FAM"
awk -F $'\t' -v OFS="\t" '{print $1,$2,"0","0","0","-9"}' "${NAME}.lgen" | sort -u > "${NAME}.fam"

# PLINK convert LGEN -> bed/bim/fam
ml plink/1.90b3.39
plink -cow --lgen "${NAME}.lgen" --fam "${NAME}.fam" --map "${NAME}.map" --recode -out "${NAME}"

# Update IDs
cp -f "$UPDATE_IDS" .
dos2unix "$(basename "$UPDATE_IDS")" || true

plink -cow -file "${NAME}" --update-ids "$(basename "$UPDATE_IDS")" --make-bed -out "UFID_${NAME}"

# Update chr / bp / alleles
cp -f "$CHR_UPD" "$BP_UPD" "$INDELS" .
dos2unix "$(basename "$CHR_UPD")" || true
dos2unix "$(basename "$BP_UPD")" || true
dos2unix "$(basename "$INDELS")" || true

plink -cow -bfile "UFID_${NAME}" --update-chr "$(basename "$CHR_UPD")" --make-bed --out mydata2
plink -cow -bfile mydata2 --update-map "$(basename "$BP_UPD")" --make-bed --out mydata3
plink -cow -bfile mydata3 --update-alleles "$(basename "$INDELS")" --recode --make-bed --out "UFID_${NAME}_Illumina.ID_ARS"

# Copy final outputs up one level
cd "$out_dir"
cp -f "./PLINK.Prep/UFID_${NAME}_Illumina.ID_ARS."* .
rm -rf ./PLINK.Prep
echo "Done."
