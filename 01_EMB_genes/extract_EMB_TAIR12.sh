#!/bin/bash

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# --- Paths and settings ---
COL=17          # Column number with gene names
FASTA_IN="TAIR12_new.faa"
FASTA_OUT="EMB_all_iso.fasta"
FASTA_LONGEST_ISO="EMB_longest_iso.fasta"
GENELIST="EMB_loci.txt"
TEMPLIST="list.txt"

# --- Step 1: Extract gene names from CSV ---
<<Comment #This only works when the linebreaks are exported as expected. In this case, proceed after the comment block
CSV= "EMB_genes.csv"
# tail -n +3 skips the first two rows (header + empty or meta row)
tail -n +3 "$CSV" \
    | cut -d',' -f"$COL" \
    | tr -d '\r' \
    | tr -d '"' \
    | grep -v '^$' > "$GENELIST" #Non-empty lines get saved to the intermediate txt file

#Idiot way to filter, if linebreaks dont match up:
#Everything gets flattened, then every 30th field after the first 17th field gets extracted
tr -d '\r\n' < "$CSV" \
    | tr ',' '\n' \
    | tr -d '"' \
    | awk 'NR == 17 || (NR > 17 && (NR - 17) % 30 == 0)' \
    | tail -n +3 \
    | grep -v '^$' > "$GENELIST"
Comment

tr -d '\r' < "$GENELIST" \
    | grep -v '^$' > "$TEMPLIST"


#echo "Genes found:$(cat < "$TEMPLIST")"
echo "Found $(wc -l < "$TEMPLIST") EMB-genes" #should be ~510 gene names

# --- Step 2: Extract matching FASTA entries ---
#NR==FNR means that this action surrounded by {} is only performed on the first file, the genelist
#/^>/ searches for fasta header
awk '
    NR == FNR {
        genes[tolower($1)] = 1 #stores genes in an associative array, dictionary. Values dont matter here
        next
    }
    /^>/ {
        match_found = 0 #resets the flag only after finding a new header
        tolower_header = tolower($0)
        for (gene in genes) {
            if (index(tolower_header, gene) > 0) {
                match_found = 1  #searches for the gene name in the current line ($0). If found, returns exit code 1
                break
            }
        }
    }
    match_found { print }
' "$GENELIST" "$FASTA_IN" > "$FASTA_OUT" #first two files are input, output is after the >

NR_EMB_GENES=$(wc -l< $TEMPLIST)
NR_FASTA_OUT=$(grep '>' < $FASTA_OUT | wc -l)

echo "Done → $(grep -c '^>' "$FASTA_OUT") matching sequences (EMB-genes and isoforms) written to $FASTA_OUT"
echo "Average number of Isoforms per EMB-gene:"
echo "scale=2 ; $NR_FASTA_OUT / $NR_EMB_GENES" | bc


# --- Step 3: Keep only longest isoform per locus ---
#awk automatically initializes new variables, they dont have to be declared seperately
awk '
    # Parse each entry, extract locus ID, track longest sequence per locus
    /^>/ {
        # Save previous entry if we have one
	# Start of LOOP
        if (header != "") {
            if (length(seq) > best_len[locus]) {
                best_len[locus]    = length(seq)
                best_header[locus] = header
                best_seq[locus]    = seq
            }
        }
        header = $0 #current line
        seq = ""
        match($0, /AT[0-9]G[0-9]{5}/, arr) #REgex searches for gene locus
        locus = toupper(arr[0])
	next
    }
    # Accumulate sequence lines, everything after current fasta header (which is skipped by the "next" statement two lines above this one) and before the new Fasta header
    { seq = seq $0 "\n" }
    # End of LOOP
    # After last entry
    END {
        if (header != "") {
            if (length(seq) > best_len[locus]) {
                best_len[locus]    = length(seq)
                best_header[locus] = header
                best_seq[locus]    = seq
            }
        } #print longest isoforms
        for (l in best_header) {
            print best_header[l]
            printf "%s", best_seq[l]
        }
    }
' "$FASTA_OUT" > "$FASTA_LONGEST_ISO"

echo "Done → $(grep -c '^>' "$FASTA_LONGEST_ISO") longest isoforms written to $FASTA_LONGEST_ISO" #Should match number of EMB-genes

rm $TEMPLIST
