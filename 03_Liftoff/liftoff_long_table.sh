#!/bin/bash
# =============================================================================
# liftoff_long_table.sh
#
# Generates a long-format table for every unique combination of:
#   gene_name x accession x copy_number
# with identity, coverage, start, and end coordinates.
#
# Output columns (TSV):
#   gene | accession | copy_number | identity | coverage | start | end
# =============================================================================

set -euo pipefail

# ---------- defaults ---------------------------------------------------------
INDIR=/ebio/abt6_projects/imgenediv/Elias/03_Liftoff/best_assembly_per_accession_liftoff
REF_GFF=best_assembly_per_accession_liftoff/100002.gff.gz #Col0
OUT="liftoff_long.tsv"
TMPDIR=$(mktemp -d)

trap 'rm -rf "$TMPDIR"' EXIT

# ---------- timing helpers bc its cool ---------------------------------------------------
SCRIPT_START=$(date +%s)
STEP_START=$SCRIPT_START

elapsed() {
    # total time since script start
    local now delta h m s
    now=$(date +%s)
    delta=$(( now - SCRIPT_START ))
    h=$(( delta / 3600 ))
    m=$(( (delta % 3600) / 60 ))
    s=$(( delta % 60 ))
    if   (( h > 0 )); then printf "%dh %dm %ds" "$h" "$m" "$s"
    elif (( m > 0 )); then printf "%dm %ds"      "$m" "$s"
    else                   printf "%ds"           "$s"
    fi
}

step_start() {
    STEP_START=$(date +%s)
}

step_elapsed() {
    # time since last step_start call
    local now delta h m s
    now=$(date +%s)
    delta=$(( now - STEP_START ))
    h=$(( delta / 3600 ))
    m=$(( (delta % 3600) / 60 ))
    s=$(( delta % 60 ))
    if   (( h > 0 )); then printf "%dh %dm %ds" "$h" "$m" "$s"
    elif (( m > 0 )); then printf "%dm %ds"      "$m" "$s"
    else                   printf "%ds"           "$s"
    fi
}

# =============================================================================
# Step 0 – build set of valid reference gene IDs from TAIR10 GFF
# =============================================================================
step_start
echo "[1/4] Extracting reference gene IDs ..."

zcat $REF_GFF | gawk '$3 == "gene" {
    n=split($9, attrs, ";")
    for (i=1; i<=n; i++) {
        split(attrs[i], kv, "=")
        if (kv[1] == "ID") print kv[2]
    }
}' - \
| gawk -v goifile="../01_EMB_genes/EMB_loci.txt" '
    BEGIN {
        while ((getline line < goifile) > 0) goi[toupper(line)] = 1
    }
    toupper($0) in goi
' > "$TMPDIR/ref_ids.txt"

REF_COUNT=$(wc -l < "$TMPDIR/ref_ids.txt")
echo "    Found $REF_COUNT reference gene IDs.  ($(step_elapsed))"


# =============================================================================
# Step 2 – collect per-accession gene records
#
# ID resolution (mirrors the pivot script logic):
#   Primary entry:   ID is directly in the whitelist
#                    -> gene = ID, copy_number = extra_copy_number + 1
#   Duplicate entry: ID has a _N suffix not in the whitelist, but base ID is
#                    -> gene = base ID, copy_number = N + 1
#
# Attributes:
#   coverage    -> coverage column
#   sequence_ID -> identity column (liftoff's non-obvious naming)
#   $4, $5      -> start, end
# =============================================================================
step_start
echo "[2/4] Scanning liftoff GFF files ..."

find "$INDIR" -name "*.gff.gz" | sort | while read -r gff; do
    ACC=$(basename "$gff" .gff.gz)

    zcat "$gff" | gawk -F'\t'\
        -v acc="$ACC" \
        -v refids="$TMPDIR/ref_ids.txt" \
    '
    BEGIN {
        while ((getline line < refids) > 0) valid[line] = 1
    }

    /^#/ { next }
    $3 == "gene" {

        id    = ""
        extra = 0
        cov   = "NA"
        ident = "NA"
        start = $4
        end   = $5

        n = split($9, attrs, ";")
        for (i = 1; i <= n; i++) {
            split(attrs[i], kv, "=")
            if      (kv[1] == "ID")                id    = kv[2]
            else if (kv[1] == "extra_copy_number") extra = kv[2] + 0
            else if (kv[1] == "coverage")          cov   = kv[2]
            else if (kv[1] == "sequence_ID")       ident = kv[2]
        }

        if (id == "") next

        # --- primary entry: ID is directly in the whitelist ---
        if (id in valid) {
            copy_number = extra + 1
            print acc "\t" id "\t" copy_number "\t" ident "\t" cov "\t" start "\t" end
        }
        # --- duplicate entry: ID has a _N suffix, base ID is in whitelist ---
        else if (id ~ /_[0-9]+$/) {
            base = id
            sub(/_[0-9]+$/, "", base)
            if (base in valid) {
                n_copy = id
                sub(/.*_/, "", n_copy)
                copy_number = n_copy
                print acc "\t" base "\t" copy_number "\t" ident "\t" cov "\t" start "\t" end
            }
        }
    }
    '
done > "$TMPDIR/raw.tsv"

LINE_COUNT=$(wc -l < "$TMPDIR/raw.tsv")
echo "    Collected $LINE_COUNT gene records.  ($(step_elapsed))"

# =============================================================================
# Step 2 – sort and write final table with header
# =============================================================================
step_start
echo "[3/4] Sorting and writing output ..."

{
    printf "gene\taccession\tcopy_number\tidentity\tcoverage\tstart\tend\n"
    # sort by gene name, then accession, then copy_number (numeric)
    sort -t$'\t' -k1,1 -k2,2 -k3,3n "$TMPDIR/raw.tsv"
} > "$OUT"

FINAL_COUNT=$(( $(wc -l < "$OUT") - 1 ))
echo "    Wrote $FINAL_COUNT rows to: $OUT  ($(step_elapsed))"

#------filter liftoff data for the EMB-genes, which weren't detected as much in the MMseq run
echo "[4/4] Filtering rows of genes with low identity in the helixer annotation"

head -1 "$OUT" > missing_genes_liftoff.tsv

grep -Ff summary_batch_run/low_presence_hits_90_90.txt "$OUT" >> missing_genes_liftoff.tsv
FILTER_COUNT=$(( $(wc -l < missing_genes_liftoff.tsv) - 1 ))

echo "    Liftoff results were filtered after the genes with low presence: $FILTER_COUNT rows written to missing_genes_liftoff.tsv"

echo
echo "Finished script: total time: $(elapsed)"
