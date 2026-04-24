#!/bin/bash

INDIR=/ebio/abt6_projects/imgenediv/Elias/02_MMseq2/best_assembly_per_accession_liftoff
OUT=liftoff_pivot.tsv
TMPDIR=$(mktemp -d)
REF_GFF=best_assembly_per_accession_liftoff/100002.gff.gz #Col0

# Step 0: build set of valid reference gene IDs from col-0 and filter only the EMB genes
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

# Step 1: collect, keeping only IDs present in reference

find "$INDIR" -name "*.gff.gz" | while read -r gff; do
    ACC=$(basename "$gff" .gff.gz)
    zcat "$gff" | gawk -F'\t' -v acc="$ACC" -v refids="$TMPDIR/ref_ids.txt" '
        BEGIN {
            while ((getline line < refids) > 0) valid[line] = 1
        }
        $3 == "gene" {
            id=""; extra=0
            n=split($9, attrs, ";")
            for (i=1; i<=n; i++) {
                split(attrs[i], kv, "=")
                if (kv[1] == "ID")                id=kv[2]
                if (kv[1] == "extra_copy_number") extra=kv[2]
            }
            if (id == "") next

            # primary entry — use extra_copy_number field
            if (id in valid) {
                total = extra + 1
                if (total > count[id]) count[id] = total
            }
            # duplicate entry — derive total from _N suffix
            else if (id ~ /_[0-9]+$/) {
                base = id
                sub(/_[0-9]+$/, "", base)
                if (base in valid) {
                    n_copy = id
                    sub(/.*_/, "", n_copy)
                    total = n_copy + 1
                    if (total > count[base]) count[base] = total
                }
            }
        }
        END {
            for (gene in count) print acc"\t"gene"\t"count[gene]
        }
    '
done > "$TMPDIR/raw.tsv"

# Step 2: pivot the temp file
gawk -F'\t' '
    {
        acc=$1; gene=$2; extra=$3
        genes[gene] = 0
        accs[acc] = 0
        val[gene][acc] = extra
    }
    END {
        n=asorti(accs, acc_sorted)
        m=asorti(genes, gene_sorted)

        printf "gene"
        for (i=1; i<=n; i++) printf "\t%s", acc_sorted[i]
        print ""

        for (j=1; j<=m; j++) {
            printf "%s", gene_sorted[j]
            for (i=1; i<=n; i++) {
                v = (acc_sorted[i] in val[gene_sorted[j]]) ? val[gene_sorted[j]][acc_sorted[i]] : "NA" #if-else shortcut (ternary statement)
                printf "\t%s", v
            }
            print ""
        }
    }
' "$TMPDIR/raw.tsv" > "$OUT"

rm -rf "$TMPDIR"

echo "liftoff count table was created: ${OUT}"
