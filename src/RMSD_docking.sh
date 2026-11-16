#!/bin/bash

PREFIX=$1   # ex: NLFR
ANALYSIS="${2:-}" # ex: WT_, or could be nothing ; THE UNDERSCORE IS IMPORTANT
LOG_DIR="/net/well/pool/projects2/s253730/projects/pred_TCR_fingerprints/logs"
mkdir -p "$LOG_DIR"

# Log file on which we will print the output of the 2 first python scripts
PREPROC_LOG="${LOG_DIR}/${ANALYSIS}${PREFIX}_dock_preproc.log"
BASE_DATA="/net/well/pool/projects2/s253730/projects/pred_TCR_fingerprints/data/${PREFIX}_analysis/af3_output_base3"

source /net/well/pool/projects2/s253730/utils/Miniconda/etc/profile.d/conda.sh
conda activate special_course

# Additionnal analysis step: compute the RMSD distance matrix between every structure, based on the docking output

python /home/projects2/pbacu/repositories/TCRdock/compute_docking_rmsds.py \
   --docking_geometries_tsvfile "${BASE_DATA}/${ANALYSIS}docking_output/parse_output.tsv" \
   --outfile "${BASE_DATA}/${ANALYSIS}docking_output/docking_rmsds.tsv" \
   >> "${LOG_DIR}/${ANALYSIS}${PREFIX}_dock_RMSD.logs" 2>&1
