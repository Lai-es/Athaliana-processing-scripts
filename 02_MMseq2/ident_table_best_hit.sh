# this step only changes the EMB gene name and does not filter for an ID threshold
 rm batch_run_mmseqs/results_longest_isoform/*_filtered.m8 || true

 for f in batch_run_mmseqs/results_longest_isoform/*.m8; do
     awk 'NR==FNR { map[$1]=$2; next } $1 in map { $1=map[$1] } 1' OFS='\t' id_mapping.tsv "$f" \
       | sort -k1,1V \
       > "${f%.m8}_filtered.m8"
 done

 echo "EMB fasta headers were renamed"


#This step pastes the best fident hit for each gene in every accession in the count table
awk 'FNR==1 {
        split(FILENAME, a, ".x.");
        sub("_filtered.m8", "", a[2]);
        accession=a[2]
        if (!(accession in col_order)) { col_order[++n_cols]=accession }
    }
    {
        if ($3 > best[$1][accession]) {
            best[$1][accession] = $3
        }
        genes[$1]=1
    }
    END {
        # print header to stderr
        printf "locus" > "/dev/stderr"
        for (i=1; i<=n_cols; i++) printf "\t%s", col_order[i] > "/dev/stderr"
        print "" > "/dev/stderr"
        # print rows to stdout
        for (gene in genes) {
            printf "%s", gene
            for (i=1; i<=n_cols; i++) printf "\t%s", (best[gene][col_order[i]] ? best[gene][col_order[i]] : 0)
            print ""
        }
    }' batch_run_mmseqs/results_longest_isoform/*_filtered.m8 2>header.tmp | sort -k1,1V > data.tmp

cat header.tmp data.tmp > summary_batch_run/count_table_best_hit.tsv
rm header.tmp data.tmp

echo "count table for the best hits is done"
