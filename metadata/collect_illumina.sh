#!/usr/bin/env bash

paths=(
"/blue/mateescu/raluca/UF Genotypes Feb2024/"
)

OUTPUT="/blue/mateescu/gzayas97/UF-LGP/metadata/finalreports.csv"

echo "full_path" > "$OUTPUT"

tmpfile=$(mktemp)

for INPUT in "${paths[@]}"; do
    if [ -d "$INPUT" ]; then
        find "$INPUT" -type f -name "*FinalReport*.txt" >> "$tmpfile"
    elif [ -f "$INPUT" ]; then
        BASENAME=$(basename "$INPUT")
        if [[ "$BASENAME" == *FinalReport*.txt ]]; then
            echo "$INPUT" >> "$tmpfile"
        fi
    else
        echo "Warning: $INPUT does not exist" >&2
    fi
done

# Only resolve paths if we actually found files
if [ -s "$tmpfile" ]; then
    sort -u "$tmpfile" >> "$OUTPUT"
else
    echo "No FinalReport files found." >&2
fi

rm -f "$tmpfile"

echo "Saved paths to $OUTPUT"