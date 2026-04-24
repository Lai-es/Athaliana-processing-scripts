#!/bin/bash
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N mmseqs
# one job uses 4 threads
#$ -pe parallel 4
# memory PER THREAD
#$ -l h_vmem=8G
# walltime per job
#$ -l h_rt=1:00:00

# ARRAY: set the range to your number of files; throttle with -tc
# e.g., 990 files total, at most 20 concurrent jobs:
#$ -t 1-460
#$ -tc 20


set -euo pipefail
set -x

# --- user params ---
CONDA_ENV=/tmp/global2/eriedling/conda/envs/mmseqs
PAIR_LIST=/ebio/abt6_projects/imgenediv/Elias/02_MMseq2/pairs_batch_run.tsv
OUT_DIR=/ebio/abt6_projects/imgenediv/Elias/02_MMseq2/batch_run_mmseqs/results_longest_isoform
TMP_DIR=/ebio/abt6_projects/imgenediv/Elias/02_MMseq2/batch_run_mmseqs/mmseq2_tmp

#make pair file
find ~+ -type f -printf "/ebio/abt6_projects/imgenediv/Elias/01_EMB_genes/EMB_longest_iso.fasta\t%p\n" > ../pairs_batch_run.tsv

mkdir -p "${OUT_DIR}"
mkdir -p "${TMP_DIR}"

TASK_ID="${SGE_TASK_ID}"

pair_line=$(sed -n "${TASK_ID}p" "$PAIR_LIST") \
    || { echo "No pair for task $TASK_ID"; exit 1; }

ref_file=$(cut -f1 <<< "$pair_line")
qry_file=$(cut -f2 <<< "$pair_line")

[[ -f "$ref_file" ]] || { echo "Missing ref file: $ref_file"; exit 1; }
[[ -f "$qry_file" ]] || { echo "Missing query file: $qry_file"; exit 1; }
#mapfile -t protein_files < <(find "$IN_PRO_DIR" -type f)
#N=${#protein_files[@]}

#task_num=$((SGE_TASK_ID - 1))
#ref_file_idx=$((task_num / N))
#query_file_idx=$((task_num % N))

#ref_file="${protein_files[$ref_file_idx]}"
#qry_file="${protein_files[$query_file_idx]}"

# set separate working directory for each job
# avoid disk I/O error
JOBTAG="${JOB_ID:-manual}.${SGE_TASK_ID:-0}"
WRK_DIR=$(mktemp -d -p "$TMP_DIR" "mmseqs_pandagma_${JOBTAG}_XXXX")
MMSEQ_TEMP="${WRK_DIR}/tmp"

mkdir -p "$MMSEQ_TEMP"

# set tmpdir
export TMPDIR="${WRK_DIR}"
export SQLITE_TMPDIR="${WRK_DIR}"
ulimit -n 4096 || true

# clean up temp folders
cleanup() { rm -rf "${WRK_DIR}" || true;
	   # rm -rf "${TMP_DIR}" || true;
		 }
trap cleanup EXIT INT TERM QUIT

# pick the input for this task
#INPUT=$(sed -n "${SGE_TASK_ID}p" "${IN_LIST}")
#[ -z "${INPUT}" ] && { echo "No input for task ${SGE_TASK_ID}"; exit 1; }

# derive unique input and output names
ref_bn=$(basename "${ref_file}")
qry_bn=$(basename "${qry_file}")
ref_sp="${ref_bn%%.*}"
qry_sp="${qry_bn%%.*}"
OUT_PATH="${OUT_DIR}/${ref_sp}.x.${qry_sp}.m8"

# align threads with allocated slots
# when NSLOTS (-pe) not given, default value will be taken (it is 8 here)
THREADS="${NSLOTS:-8}"

# run
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"
cd "$WRK_DIR"
mmseqs easy-search \
    "$ref_file" "$qry_file" "$OUT_PATH" "$MMSEQ_TEMP" \
    --search-type 1 --cov-mode 0 \
    --threads 4 --split-memory-limit 6G \
    --format-output "query,target,fident,qcov,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln" 1>/dev/null
