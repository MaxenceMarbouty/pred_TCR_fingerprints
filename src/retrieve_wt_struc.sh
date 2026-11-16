#!/bin/bash

PREFIX=$1   # ex: NLFR

AF3_BASE_DIR="/net/well/pool/projects2/s253730/projects/pred_TCR_fingerprints/results/${PREFIX}_analysis/af3_output/structInference_${PREFIX}/json_unpairedMSA_onqueryTemplate/${PREFIX}-P0_WT"
DEST_DIR="../data/${PREFIX}_analysis/af3_output_base3/WT_structures"


mkdir -p "$DEST_DIR"

count_ok=0
count_missing=0

for seed_nb in {1..10}; do

    for sample_nb in {0..4}; do

        # Reconstructing name of cif model
        cif_filename="${PREFIX}-P0_WT_seed-${seed_nb}_sample-${sample_nb}_model.cif"
        src_file="${AF3_BASE_DIR}/seed-${seed_nb}_sample-${sample_nb}/${cif_filename}"

        # Verifying the file exists before copying it 
        if [ -f "$src_file" ]; then
            cp "$src_file" "${DEST_DIR}/seed-${seed_nb}_sample-${sample_nb}.cif"
            
            ((count_ok++))
        else
            echo "Missing file: $src_file"
            ((count_missing++))
        fi
    done
done

echo "Structures successfully copied: $count_ok"
echo "Missing files : $count_missing"
echo "Expected total : 50"