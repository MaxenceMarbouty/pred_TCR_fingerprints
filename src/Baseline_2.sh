#!/bin/bash

PREFIX=$1   # ex: NLFR
HLA=$2      # ex: HLA-B08:01
EXP=$3      # ex: TCRneo

tail -n +2 ../data/${PREFIX}_analysis/${EXP}_duplex_fingerprint_data_WT_TCR.tsv | cut -f5 | sort -u > ../data/${PREFIX}_analysis/peptides.txt

/home/ctools/netMHCpan-4.2/netMHCpan -p ../data/${PREFIX}_analysis/peptides.txt -a ${HLA} > ../results/netmhcpan_${PREFIX}_analysis.txt

/home/ctools/netMHCpan-4.2/netMHCpan -p ../data/${PREFIX}_analysis/peptides.txt -a ${HLA} -xls -xlsfile ../results/${PREFIX}_analysis/netmhcpan_${PREFIX}_analysis_clean.tsv
