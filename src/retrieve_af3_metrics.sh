#!/bin/bash

PREFIX=$1   # ex: NLFR


MAPPING_CSV="../data/${PREFIX}_analysis/af3_output_base3/selected_structures.csv"
AF3_BASE_DIR="/net/well/pool/projects2/s253730/projects/pred_TCR_fingerprints/results/${PREFIX}_analysis/af3_output/structInference_${PREFIX}/json_unpairedMSA_onqueryTemplate"
DEST_DIR="../data/${PREFIX}_analysis/af3_output_base3/metrics"


mkdir -p "$DEST_DIR"


count_ok=0
count_missing=0

tail -n +2 "$MAPPING_CSV" | while IFS=',' read -r pdb_id seed_sample rest_of_line; do
    
    pdb_id=$(echo "$pdb_id" | tr -d '\r')
    seed_sample=$(echo "$seed_sample" | tr -d '\r')

    # Reconstructing name of cif model
    json_filename="${pdb_id}_${seed_sample}_confidences.json"
    src_file="${AF3_BASE_DIR}/${pdb_id}/${seed_sample}/${json_filename}"

    # Verifying the file exists before copying it 
    if [ -f "$src_file" ]; then
        cp "$src_file" "${DEST_DIR}/${pdb_id}_confidences.json"
        
        ((count_ok++))
    else
        echo "Missing file: $src_file"
        ((count_missing++))
    fi

done < <(tail -n +2 "$MAPPING_CSV")

echo "Structures successfully copied: $count_ok"
echo "Missing files : $count_missing"
echo "Expected total : $(tail -n +2 "$MAPPING_CSV" | wc -l)"