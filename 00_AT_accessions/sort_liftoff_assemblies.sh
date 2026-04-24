#!/bin/bash

# folders
CSV="metadata_athaliana.csv"
FASTA_DIR="./liftoff_raw"
OUT_DIR="./best_assembly_per_accession_liftoff"

mkdir -p "${OUT_DIR}"

# Column numbering to cut later
COL_ASSEMBLY_ID=1
COL_ACCESSION=2
COL_BEST=19

# Skip header, extract only relevant columns with CUT
tail -n +2 "$CSV" \
    | cut -d',' -f"${COL_ASSEMBLY_ID},${COL_ACCESSION},${COL_BEST}" \
    | while IFS=',' read -r assembly_id accession best; do
    # Whitespace cleanup
    assembly_id=$(echo "$assembly_id" | tr -d '\r' | xargs)
    accession=$(echo "$accession" | tr -d '\r' | xargs)
    best=$(echo "$best" | tr -d '\r' | xargs)
    if [[ "$best" == "yes" ]]; then
        src="$FASTA_DIR/${assembly_id}.TAIR12.genes.gff.gz"
        dst="$OUT_DIR/${assembly_id}.gff.gz" #change to ${accession} if file should be named like Col-0.pro.fa.gz etc...
        if [[ -f "$src" ]]; then
            cp "$src" "$dst"
            echo " ^|^s $assembly_id  ^f^r $accession.fasta"
        else
            echo " ^z  Datei nicht gefunden: $src"
        fi
    fi
done
