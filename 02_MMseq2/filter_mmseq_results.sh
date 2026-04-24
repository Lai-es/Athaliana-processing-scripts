#!/bin/bash

#----------------------User params-------------------
OUTDIR=./summary_batch_run #output directory
ID=90 #identity threshhold in percent
COV=90 #coverage treshhold in percent

IDp=0.9
COVp=0.9
#--------------------------------------------------------------------------Step 1------------------------------------
# This step extracts the locus id from the EMB-fasta file and replaces it in the MMseq2-output files

#Extracts fasta headers and writes a two-columns tsv. first column: start of the fasta header. second column: gene locus id
grep "^>" ../01_EMB_genes/EMB_all_iso.fasta | sed 's/^>\([^ ]*\).*\[locus_tag=TAIR12_TAIR12_\([^]]*\)\].*/\1\t\2/; s/^lcl|//' > id_mapping.tsv

echo "Mapping file was created"

#generalization for multiple files. For every file, if the start of the fasta header is found in the id_mapping tsv, substitute it with the gene locus id

rm batch_run_mmseqs/results_longest_isoform/*_filtered.m8 || true #cleanup so it doesnt run for the already filtered versions

for f in batch_run_mmseqs/results_longest_isoform/*.m8; do
    awk 'NR==FNR { map[$1]=$2; next } $1 in map { $1=map[$1] } 1' OFS='\t' id_mapping.tsv "$f" \
      | awk -v ID=$IDp '$3 > ID' \
      | awk -v COV=$COVp '$4 > COV' \
      | sort -k1,1V \
      > "${f%.m8}_filtered.m8"
done
#two sequential filtering steps to look first if there is a match with ID $3 > 90%, which also have coverage $4 > 90%

echo "MMseq2 results were filtered for ID > ${ID}% and coverage > ${COV}%, EMB fasta headers were renamed"

#------------------------------------------------------------------Step 2-----------------------------------
# In this step, a table is created for every run to investigate how often each emb gene is present in each accession with ID>0.90

awk 'FNR==1 {
        # extract accession: take part after ".x.", strip "_filtered"
        split(FILENAME, a, ".x.");
        sub("_filtered.m8", "", a[2]);
        accession=a[2] #eg Col-0
        if (!(accession in col_order)) { col_order[++n_cols]=accession }
    }
    {
        count[$1][accession]++ #$1 is the gene locus, $3 is fidemt, the fraction of identical sites
        genes[$1]=1
    }
    END {
        # print header
        printf "EMB-locus" > "/dev/stderr"
        for (i=1; i<=n_cols; i++) printf "\t%s", col_order[i] > "/dev/stderr"
        print "" > "/dev/stderr"
        # print rows
        for (gene in genes) {
            printf "%s", gene
            for (i=1; i<=n_cols; i++) printf "\t%s", (count[gene][col_order[i]] ? count[gene][col_order[i]] : 0)
            print ""
        }
    }' batch_run_mmseqs/results_longest_isoform/*_filtered.m8 2>header.tmp | sort -k1,1V > table.tmp

cat header.tmp table.tmp > ${OUTDIR}/count_table_${ID}_${COV}.tsv
rm header.tmp table.tmp

echo "A count table was made in which the appearance with ID>90% of EMB genes was counted for each accession"

#-----------------------------------------------------------Step 3--------------------------------
#Generate a histogram from the data

awk 'NR>1 { for (i=2; i<=NF; i++) freq[$i]++ }
     END { for (val in freq) print val"\t"freq[val] }' summary_batch_run/count_table_${ID}_${COV}.tsv | sort -k1,1n > ${OUTDIR}/histogram_${ID}_${COV}.tsv


echo "Histogram for the mmseq2 results:"

awk -v width="$(tput cols)" '
    NR==FNR { if ($2>max) max=$2; next }
    {
        bar_width = width - 20
        len = int($2 / max * bar_width)
        bar = ""
        for (i=0; i<len; i++) bar = bar "█"
        printf "%5s | %-*s %s\n", $1, bar_width, bar, $2
    }' ${OUTDIR}/histogram_${ID}_${COV}.tsv  ${OUTDIR}/histogram_${ID}_${COV}.tsv


#---------------------------------------------------Step 4------------------------------------
#generate presence fraction for each gene

awk 'NR==1 {
        n_acc = NF - 1
        print "locus\tpresence_fraction\tlow_presence?"
        next
    }
    {
        present = 0
        for (i=2; i<=NF; i++) if ($i > 0) present++
	frac=present/n_acc
	if (frac < 0.95)
            printf "%s\t%.4f\tyes\n", $1, frac
	else
            printf "%s\t%.4f\tno\n", $1, frac
    }' ${OUTDIR}/count_table_${ID}_${COV}.tsv >  ${OUTDIR}/presence_fractions_${ID}_${COV}.tsv


#--------------------------------------------Step 5----------------------
#Extract genes with lower than 95% presence in the presence_fractions.tsv

awk -F'\t' '$3=="yes" {print $1}' ${OUTDIR}/presence_fractions_${ID}_${COV}.tsv > ${OUTDIR}/low_presence_hits_${ID}_${COV}.txt

#---------------------Step 6-----------------------
#make a per-gene histogram table

awk '
    NR==1 { n_acc=NF-1; next }
    {
        gene = $1
        delete freq
        for (i=2; i<=NF; i++) freq[$i]++
        for (val in freq) {
            counts[gene][val] = freq[val]
            seen_vals[val] = 1
        }
        genes[++n_genes] = gene
    }
    END {
        # collect and sort observed values
        n_vals = asorti(seen_vals, sorted_vals, "@ind_num_asc")

        # header
        printf "locus"
        for (v=1; v<=n_vals; v++) printf "\t%s", sorted_vals[v]
        print ""

        # rows
        for (g=1; g<=n_genes; g++) {
            printf "%s", genes[g]
            for (v=1; v<=n_vals; v++)
                printf "\t%s", (counts[genes[g]][sorted_vals[v]] ? counts[genes[g]][sorted_vals[v]] : 0)
            print ""
        }
    }
' ${OUTDIR}/count_table_${ID}_${COV}.tsv > ${OUTDIR}/gene_copy_freq_${ID}_${COV}.tsv

echo "Per-gene copy number frequency table written to ${OUTDIR}/gene_copy_freq_${ID}_${COV}.tsv"
