#!/bin/bash
#SBATCH --time=48:00:00

proj_env=$1
NAME=$2
Illumina_meta_csv=$3

source "$proj_env"

out_dir="${my_results}/${NAME}"
prep_dir="${out_dir}/Illumina"
mkdir -p "$out_dir" "$prep_dir"

echo "Making directories for new genotypes:"
echo "  out_dir:  $out_dir"
echo "  prep_dir: $prep_dir"
echo "Copying genotypes listed in meta CSV: $Illumina_meta_csv"

# --- Copy files listed in the meta CSV ---
# Expect a header column named: full_path
# Example:
# full_path
# /blue/.../Univ_of_Florida_..._FinalReport.txt

if [[ ! -f "$Illumina_meta_csv" ]]; then
  echo "ERROR: meta CSV not found: $Illumina_meta_csv" >&2
  exit 1
fi

# Read the header, find which column is "full_path"
header="$(head -n 1 "$Illumina_meta_csv")"
IFS=',' read -r -a cols <<< "$header"

full_path_col=-1
for i in "${!cols[@]}"; do
  colname="${cols[$i]}"
  colname="${colname//$'\r'/}"      # strip CR if file is Windows-formatted
  if [[ "$colname" == "full_path" ]]; then
    full_path_col=$((i+1))          # awk is 1-based
    break
  fi
done

if [[ "$full_path_col" -lt 1 ]]; then
  echo "ERROR: Could not find column named 'full_path' in $Illumina_meta_csv" >&2
  echo "Header was: $header" >&2
  exit 1
fi

# Extract full_path values (skip header), copy each file
# - Handles quoted CSV poorly (simple CSV only). If your paths can contain commas/quotes, say so.
awk -F',' -v c="$full_path_col" 'NR>1 {print $c}' "$Illumina_meta_csv" \
  | sed 's/^ *//; s/ *$//; s/\r$//' \
  | grep -v '^$' \
  | while IFS= read -r src; do
      if [[ ! -f "$src" ]]; then
        echo "WARNING: missing file (skipping): $src" >&2
        continue
      fi
      echo "  cp: $src"
      cp -f "$src" "$prep_dir/"
    done

echo "Done copying. Files now in: $prep_dir"
echo

# --- Quick checks (run on the copied files) ---
echo "Preview headers:"
head -20 "$prep_dir"/*.txt || true
echo

echo "Num Samples lines found:"
grep -H "Num Samples" "$prep_dir"/*.txt || true
echo

echo "Column count sanity check (first 20 lines per file):"
shopt -s nullglob
for file in "$prep_dir"/*.txt; do
  echo "FILE: $(basename "$file")"
  head -20 "$file" | awk '{print NF}' | sort -nu | head -n 1
  head -20 "$file" | awk '{print NF}' | sort -nu | tail -n 2
  echo
done
