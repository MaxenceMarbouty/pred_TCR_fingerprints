#!/bin/bash

#SBATCH --job-name=TCRdock
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=150:00:00

#Strict execution mode
set -Eeuo pipefail
source /net/well/pool/projects2/s253730/utils/Miniconda/etc/profile.d/conda.sh 
conda activate special_course

PREFIX=$1
ANALYSIS="${2:-}"

base_folder=/net/well/pool/projects2/s253730/projects/pred_TCR_fingerprints/data/${PREFIX}_analysis/af3_output_base3/${ANALYSIS}structures_pdb_cleaned 
output_folder=/net/well/pool/projects2/s253730/projects/pred_TCR_fingerprints/data/${PREFIX}_analysis/af3_output_base3/${ANALYSIS}docking_output
mkdir -p $output_folder

pdb_files=$(ls $base_folder/*.pdb | tr '\n' ' ')
#exclude 7rm4.trunc.fit.pdb from the list of pdb files to parse, as it is a truncated version of 7rm4.pdb
# pdb_files=$(echo $pdb_files | tr ' ' '\n' | grep -v "7rm4.trunc.fit.pdb" | tr '\n' ' ')
python /home/projects2/pbacu/repositories/TCRdock/parse_tcr_pmhc_pdbfile.py --pdbfiles $pdb_files \
    --organism human --mhc_class 1 --out_tsvfile $output_folder/parse_output.tsv