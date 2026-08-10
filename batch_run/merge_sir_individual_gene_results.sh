#!/bin/bash

clearance_gene_dir="../../derivatives/yohan_ge_filt/mr10vv0.2/"
clearance_list_path="${clearance_gene_dir}gene_names_mr10vv0.2.csv"

# Load full gene list (assuming one gene per line)
mapfile -t all_genes < "$clearance_list_path"

find ../sir_result_csvs/clearance/ -mindepth 1 -type d | while read -r subdir; do

    # Remove existing output files if present
    [[ -f "$subdir/final_output.csv" ]] && rm "$subdir/final_output.csv" && echo "Removed existing final_output.csv"
    [[ -f "$subdir/missing_genes.csv" ]] && rm "$subdir/missing_genes.csv" && echo "Removed existing missing_genes.csv"
    [[ -f "$subdir/top40_genes.csv" ]] && rm "$subdir/top40_genes.csv" && echo "Removed existing top40_genes.csv"

    files_present=$(find "$subdir" -maxdepth 1 -type f)
    [[ -z "$files_present" ]] && echo "Skipping (no relevant files): $subdir" && continue

    echo "=== Processing: $subdir ==="

    shopt -s nullglob
    files=("$subdir"/*.csv)

    nonempty=()

    # Collect non-empty csv files (excluding final_output.csv and missing_genes.csv)
    for f in "${files[@]}"; do
        [[ "$f" == "$subdir/final_output.csv" ]] && continue
        [[ "$f" == "$subdir/missing_genes.csv" ]] && continue
        [[ "$f" == "$subdir/top40_genes.csv" ]] && continue
        [[ -s "$f" ]] && nonempty+=("$f")
    done

    # Extract gene names from non-empty file basenames, sort them
    found_genes=$(printf '%s\n' "${nonempty[@]}" | xargs -I{} basename {} .csv | sort)

    # Compare sorted lists: comm -23 shows lines only in all_genes (i.e. missing from found)
    missing_genes=$(comm -23 <(printf '%s\n' "${all_genes[@]}" | sort) <(echo "$found_genes"))

    # Merge non-empty files, keeping best row per gene (highest col 2)
    if [[ ${#nonempty[@]} -gt 0 ]]; then
        echo "Merging ${#nonempty[@]} files in: $subdir"
        cat "${nonempty[@]}" | awk -F',' 'NF > 1 && ($2+0) > seen[$1] { seen[$1]=$2; best[$1]=$0 } END { for (g in best) print best[g] }' >> "$subdir/final_output.csv"

        # Create top40_genes.csv by sorting final_output.csv by col 2 descending, taking top 40
        { echo "gene"; sort -t',' -k2,2rn "$subdir/final_output.csv" | head -n40 | awk -F',' '{print $1}'; } > "$subdir/top40_genes.csv"
    else
        echo "Skipping (no non-empty csv files): $subdir"
    fi

    # Report missing genes
    if [[ -n "$missing_genes" ]]; then
        echo "Missing genes in $subdir:"
        echo "$missing_genes"
        { echo "gene"; echo "$missing_genes"; } > "$subdir/missing_genes.csv"
    else
        echo "All genes accounted for in $subdir"
    fi

done