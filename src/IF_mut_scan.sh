#!/bin/bash
#SBATCH --job-name=esm_if_scan
#SBATCH --output=../logs/SLLQ_analysis/esm_if_%j.log
#SBATCH --error=../logs/SLLQ_analysis/esm_if_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1


source /net/well/pool/projects2/s253730/utils/Miniconda/etc/profile.d/conda.sh 
conda activate inverse_folding

# Don't forget to change the .log and .err path as well 
PEP_FOLDER="SLLQ_analysis"
PEP_SHORT="SLLQ"
PDB_CLEANED="../data/${PEP_FOLDER}/af3_output_base3/structures_pdb_cleaned/${PEP_SHORT}-P0_WT.pdb"
FASTA_IN="../data/${PEP_FOLDER}/mutated_peptides.fasta"
CSV_OUT="../results/${PEP_FOLDER}/esm_if_scores.csv"

# Exécution du script officiel
python ./score_log_likelihoods.py \
    "${PDB_CLEANED}" \
    "${FASTA_IN}" \
    --chain C \
    --multichain-backbone \
    --outpath "${CSV_OUT}"
