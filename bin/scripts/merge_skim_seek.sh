#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1

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

proj_env=$1
NAME=$2
meta_csv=$3
# -------------------- Confidence toggles --------------------
# 1 = keep LOWCONF and low-GP genotypes as called
# 0 = set them to missing GT=./.
KEEP_LOWCONF="${$4:-0}"
# Used only when KEEP_LOWCONF=0
GP_MIN="${$5:-0.90}"
# Final variant IDs
SNP_ID_FMT="${SNP_ID_FMT:-%CHROM:%POS:%REF:%ALT}"

source $proj_env
results="${my_results}"
out_dir="${results}/${NAME}"
out_name="merged_output.vcf.gz"

NTHREADS="${SLURM_CPUS_PER_TASK:-12}"

# -------------------- Confidence toggles --------------------
# 1 = keep LOWCONF and low-GP genotypes as called
# 0 = set them to missing GT=./.
KEEP_LOWCONF="${KEEP_LOWCONF:-0}"
# Used only when KEEP_LOWCONF=0
GP_MIN="${GP_MIN:-0.90}"

# Final variant IDs
SNP_ID_FMT="${SNP_ID_FMT:-%CHROM:%POS:%REF:%ALT}"

mkdir -p "${out_dir}"
cd "${out_dir}"

# -------------------- Modules --------------------
ml bcftools
ml htslib
ml samtools

# -------------------- Helpers --------------------
die() { echo "ERROR: $*" >&2; exit 1; }
msg() { echo "INFO: $*"; }
have_parallel() { command -v parallel >/dev/null 2>&1; }

run_parallel_template() {
  # usage: run_parallel_template 'command with {} placeholder' file_list
  local tmpl="$1"
  local flist="$2"
  if have_parallel; then
    parallel -j "${NTHREADS}" "${tmpl}" :::: "${flist}"
  else
    xargs -a "${flist}" -P "${NTHREADS}" -I{} bash -c "${tmpl}"
  fi
}

need_setGT_if_enabled() {
  if [[ "${KEEP_LOWCONF}" -eq 0 ]]; then
    if ! bcftools plugin -l 2>/dev/null | grep -qx 'setGT'; then
      die "bcftools plugin 'setGT' not available in this environment. Load a bcftools build with plugins or set KEEP_LOWCONF=1."
    fi
  fi
}

# ============================================================
#  Temp workspace (prefer node-local scratch if available)
# ============================================================
scratch_root="${out_dir}"
tmp_root="${scratch_root}/skimseek_merge_${SLURM_JOB_ID:-$$}"
clean_dir="${tmp_root}/clean_vcfs"
logs_dir="${tmp_root}/logs"
merge_dir="${tmp_root}/merge_by_contig"
mkdir -p "${clean_dir}" "${logs_dir}" "${merge_dir}"

# cleanup() { rm -rf "${tmp_root}"; }
# trap cleanup EXIT

msg "Temp workspace: ${tmp_root}"
msg "KEEP_LOWCONF=${KEEP_LOWCONF}  (0=set missing, 1=keep as-is)"
msg "GP_MIN=${GP_MIN}              (used only if KEEP_LOWCONF=0)"

need_setGT_if_enabled

# ============================================================
#  1) Read metadata CSV and build TSV + path list (awk version)
# ============================================================
[[ -f "${meta_csv}" ]] || die "Metadata CSV not found: ${meta_csv}"

samples_tsv="${out_dir}/samples_2026.tsv"   # path <tab> SkimSeek_ID <tab> NewID
paths_list="${out_dir}/vcf_paths.list"      # one VCF path per line

: > "${samples_tsv}"
: > "${paths_list}"

msg "Parsing metadata (awk) -> ${samples_tsv}"

awk -v OUT_TSV="${samples_tsv}" -v OUT_LIST="${paths_list}" '
  function trim(s) { gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s); return s }
  function dequote(s) { s = trim(s); gsub(/^"/, "", s); gsub(/"$/, "", s); return s }

  BEGIN { FS=","; OFS="\t"; n_out=0; err=0 }

  NR==1 {
    for (i=1; i<=NF; i++) { h=dequote($i); col[h]=i }
    if (!("path" in col) || !("SkimSeek_ID" in col) || !("NewID" in col)) {
      print "ERROR: metadata CSV missing required columns. Need: path,SkimSeek_ID,NewID" > "/dev/stderr"
      print "ERROR: header seen: " $0 > "/dev/stderr"
      exit 2
    }
    next
  }

  NR>1 {
    path = dequote($(col["path"]))
    sid  = dequote($(col["SkimSeek_ID"]))
    nid  = dequote($(col["NewID"]))
    if (path=="" || sid=="" || nid=="") next

    if (++seen_sid[sid] > 1) { print "ERROR: Duplicate SkimSeek_ID in metadata: " sid > "/dev/stderr"; err=1 }
    if (++seen_nid[nid] > 1) { print "ERROR: Duplicate NewID in metadata: " nid > "/dev/stderr"; err=1 }

    print path, sid, nid >> OUT_TSV
    print path >> OUT_LIST
    n_out++
  }

  END {
    if (n_out==0) { print "ERROR: metadata produced 0 usable rows (need non-empty path/SkimSeek_ID/NewID)" > "/dev/stderr"; exit 3 }
    if (err) exit 4
  }
' "${meta_csv}"

[[ -s "${samples_tsv}" ]] || die "samples TSV is empty: ${samples_tsv}"
[[ -s "${paths_list}"  ]] || die "paths list is empty: ${paths_list}"
msg "Samples loaded: $(wc -l < "${samples_tsv}")"

# Check that paths exist
msg "Checking that VCF paths exist..."
missing_paths=0
while IFS= read -r p; do
  if [[ ! -f "${p}" ]]; then
    echo "WARN: Missing VCF file: ${p}"
    ((missing_paths++)) || true
  fi
done < "${paths_list}"
[[ "${missing_paths}" -eq 0 ]] || die "Missing ${missing_paths} VCFs listed in metadata."

# ============================================================
#  2) Index inputs (if needed)
# ============================================================
msg "Indexing input VCFs (if needed) with up to ${NTHREADS} workers..."
run_parallel_template '[[ -f "{}.tbi" ]] || tabix -p vcf "{}"' "${paths_list}"

# ============================================================
#  3) Safety checks: each VCF is single-sample
# ============================================================
msg "Validating each input VCF is single-sample..."
while IFS= read -r p; do
  ns=$(bcftools query -l "${p}" | wc -l)
  if [[ "${ns}" -ne 1 ]]; then
    die "${p} has ${ns} samples (expected 1)."
  fi
done < "${paths_list}"
msg "Single-sample check passed."

# ============================================================
#  4) Cleaning (parallel), and force sample name to SkimSeek_ID
#     - Keep PASS and LOWCONF records
#     - Optionally set low-confidence genotypes to missing based on FILTER/GP
#     - Strip to GT-only for merge stability
# ============================================================

clean_one_vcf() {
  local in_vcf="$1"
  local skim_id="$2"

  local out_vcf="${clean_dir}/${skim_id}.clean.vcf.gz"
  local tmp_clean="${out_vcf}.clean.tmp.gz"
  local tmp_named="${out_vcf}.named.tmp.gz"
  local log="${logs_dir}/${skim_id}.clean.log"

  if [[ -s "${out_vcf}" && -s "${out_vcf}.tbi" ]]; then
    echo "${out_vcf}"
    return 0
  fi

  # We keep both PASS and LOWCONF so we can toggle behavior later
  local keep_records_expr='(FILTER="PASS" || FILTER="LOWCONF")'

  # If we are setting to missing, define low-confidence genotypes
  # Note: MAX(FMT/GP) works because GP is a vector per sample.
  local set_missing_expr="(FILTER=\"LOWCONF\" || MAX(FMT/GP)<${GP_MIN})"

  {
    echo "INPUT:           ${in_vcf}"
    echo "SKIMSEEK_ID:     ${skim_id}"
    echo "OUTPUT:          ${out_vcf}"
    echo "KEEP_LOWCONF:    ${KEEP_LOWCONF}"
    echo "GP_MIN:          ${GP_MIN}"
    echo "KEEP_RECORDS:    ${keep_records_expr}"
    echo "SET_MISSING_IF:  ${set_missing_expr}"
    echo
    echo "Steps:"
    echo "  1) drop missing GT records                   : bcftools view -U"
    echo "  2) keep records with FILTER PASS or LOWCONF  : bcftools view -i"
    if [[ "${KEEP_LOWCONF}" -eq 0 ]]; then
      echo "  3) set low-confidence GT to missing (./.)    : bcftools +setGT"
    else
      echo "  3) keep low-confidence GT as-is              : (skip setGT)"
    fi
    echo "  4) strip INFO and keep GT only               : bcftools annotate -x ..."
    echo "  5) force sample name = SkimSeek_ID           : bcftools reheader"
  } > "${log}"

  if [[ "${KEEP_LOWCONF}" -eq 1 ]]; then
    # Keep genotypes as called; just keep PASS/LOWCONF records and strip to GT-only
    bcftools view -U -Ou "${in_vcf}" | \
      bcftools view -i "${keep_records_expr}" -Ou | \
      bcftools annotate -x INFO,FORMAT/RC,FORMAT/AC,FORMAT/GP,FORMAT/DS -Oz -o "${tmp_clean}" >> "${log}"
  else
    # Convert low-confidence genotypes to missing, then strip to GT-only
    bcftools view -U -Ou "${in_vcf}" | \
      bcftools view -i "${keep_records_expr}" -Ou | \
      bcftools +setGT -Ou -- -t q -n . -i "${set_missing_expr}" | \
      bcftools annotate -x INFO,FORMAT/RC,FORMAT/AC,FORMAT/GP,FORMAT/DS -Oz -o "${tmp_clean}" >> "${log}"
  fi

  tabix -p vcf "${tmp_clean}" 2>> "${log}"

  # Force the sample name to the SkimSeek_ID (single-sample file)
  echo "${skim_id}" > "${tmp_root}/.${skim_id}.name"
  bcftools reheader -s "${tmp_root}/.${skim_id}.name" -o "${tmp_named}" "${tmp_clean}" 2>> "${log}"
  tabix -p vcf "${tmp_named}" 2>> "${log}"

  mv -f "${tmp_named}" "${out_vcf}"
  mv -f "${tmp_named}.tbi" "${out_vcf}.tbi"

  rm -f "${tmp_clean}" "${tmp_clean}.tbi" "${tmp_root}/.${skim_id}.name"

  echo "${out_vcf}"
}

export -f clean_one_vcf
export clean_dir logs_dir tmp_root KEEP_LOWCONF GP_MIN

msg "Cleaning VCFs in parallel..."
clean_list="${out_dir}/vcf_to_merge.cleaned.list"
: > "${clean_list}"

jobs_tsv="${tmp_root}/clean_jobs.tsv"
cut -f1-2 "${samples_tsv}" > "${jobs_tsv}"

if have_parallel; then
  parallel -j "${NTHREADS}" --colsep '\t' clean_one_vcf {1} {2} :::: "${jobs_tsv}" >> "${clean_list}"
else
  xargs -a "${jobs_tsv}" -P "${NTHREADS}" -I{} bash -lc '
    IFS=$'\''\t'\'' read -r p sid <<< "{}"
    clean_one_vcf "$p" "$sid"
  ' >> "${clean_list}"
fi

[[ -s "${clean_list}" ]] || die "Cleaning produced an empty list: ${clean_list}"
msg "Cleaned VCFs ready: $(wc -l < "${clean_list}")"

# ============================================================
#  5) Merge cleaned VCFs by contig in parallel, then concat
# ============================================================
msg "Building contig list from data (includes NKLS-like contigs)..."

contigs_raw="${merge_dir}/contigs.raw.txt"
contigs_ordered="${merge_dir}/contigs.ordered.txt"

head -n 20 "${clean_list}" | while read -r v; do
  bcftools index -s "${v}" | awk 'NF{print $1}'
done | sort -u > "${contigs_raw}"

[[ -s "${contigs_raw}" ]] || die "Failed to infer contigs from cleaned VCFs."

# Order: 1..29, X, MT, then everything else
awk '
  BEGIN { OFS="\t" }
  {
    c=$1
    if (c ~ /^[0-9]+$/) {
      n=c+0
      if (n>=1 && n<=29) printf("0\t%02d\t%s\n", n, c)
      else               printf("3\t%999d\t%s\n", n, c)
    } else if (c=="X")  printf("1\t00\t%s\n", c)
    else if (c=="MT")  printf("2\t00\t%s\n", c)
    else               printf("4\t00\t%s\n", c)
  }
' "${contigs_raw}" | sort -k1,1n -k2,2n -k3,3 | awk '{print $3}' > "${contigs_ordered}"

msg "Contigs to merge: $(wc -l < "${contigs_ordered}")"

export clean_list merge_dir

merge_one_contig() {
  local c="$1"
  local out="${merge_dir}/merged.${c}.vcf.gz"
  local tmp="${out}.tmp.gz"

  if [[ -s "${out}" && -s "${out}.tbi" ]]; then
    return 0
  fi

  bcftools merge --force-samples -r "${c}" -Oz --threads 1 --file-list "${clean_list}" -o "${tmp}"
  tabix -p vcf "${tmp}"

  mv -f "${tmp}" "${out}"
  mv -f "${tmp}.tbi" "${out}.tbi"
}
export -f merge_one_contig

msg "Merging per contig in parallel..."
if have_parallel; then
  parallel -j "${NTHREADS}" merge_one_contig :::: "${contigs_ordered}"
else
  xargs -a "${contigs_ordered}" -P "${NTHREADS}" -I{} bash -lc 'merge_one_contig "$0"' {}
fi

parts_list="${merge_dir}/parts.ordered.list"
: > "${parts_list}"
while read -r c; do
  echo "${merge_dir}/merged.${c}.vcf.gz" >> "${parts_list}"
done < "${contigs_ordered}"

merged_vcf="${out_dir}/${out_name}"
msg "Concatenating contig pieces -> ${merged_vcf}"
bcftools concat -a -Oz --threads "${NTHREADS}" -f "${parts_list}" -o "${merged_vcf}"
tabix -@ "${NTHREADS}" -p vcf "${merged_vcf}"

# ============================================================
#  6) Rename samples to NewID using metadata (SkimSeek_ID -> NewID)
# ============================================================
msg "Renaming samples to NewID using metadata..."

cur_samples="${out_dir}/_current_samples.txt"
new_samples="${out_dir}/_new_samples.txt"
bcftools query -l "${merged_vcf}" > "${cur_samples}"

awk '
  NR==FNR { map[$2]=$3; next }
  {
    if ($1 in map && map[$1] != "") print map[$1];
    else                           print $1;
  }
' "${samples_tsv}" "${cur_samples}" > "${new_samples}"

[[ -s "${new_samples}" ]] || die "Failed to build NewID rename list."

renamed_vcf="${out_dir}/renamed_${out_name}"
bcftools reheader --threads "${NTHREADS}" -s "${new_samples}" -o "${renamed_vcf}" "${merged_vcf}"
tabix -@ "${NTHREADS}" -p vcf "${renamed_vcf}"

# ============================================================
#  7) Assign SNP IDs = CHR:POS:REF:ALT
# ============================================================
msg "Assigning variant IDs as ${SNP_ID_FMT} ..."

final_vcf="${out_dir}/${NAME}_${out_name}"
bcftools annotate --threads "${NTHREADS}" --set-id "${SNP_ID_FMT}" -Oz -o "${final_vcf}" "${renamed_vcf}"
tabix -@ "${NTHREADS}" -p vcf "${final_vcf}"
bcftools stats "${final_vcf}" > "${out_dir}/${NAME}_stats.txt"

msg "Final output: ${final_vcf}"
msg "DONE."
