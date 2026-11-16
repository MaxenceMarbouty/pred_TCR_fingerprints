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

# First step: Converting AF3 selected structures (one per mutation) from .cif to .pdb format

python docking_analysis/cif2pdb.py -i "${BASE_DATA}/${ANALYSIS}structures" -o "${BASE_DATA}/${ANALYSIS}structures_pdb" \
> "${PREPROC_LOG}" 2>&1

# Second step: Cleaning the resulting PdB files (removing solvent and renumbering the elements)

python docking_analysis/pdbCleaner.py "${BASE_DATA}/${ANALYSIS}structures_pdb" "${BASE_DATA}/${ANALYSIS}structures_pdb_cleaned" \
>> "${PREPROC_LOG}" 2>&1

# Last step: perform docking analysis on tailored pdb files. Requires a lot of ressources so we use the queuing system

sbatch --output="${LOG_DIR}/${ANALYSIS}${PREFIX}_docking_%j.out" \
  --error="${LOG_DIR}/${ANALYSIS}${PREFIX}_docking_%j.err" \
  docking_analysis/launchTCRdock_templates.sh "${PREFIX}" "${ANALYSIS}"

