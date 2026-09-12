
tail -n +2 ../data/TCRneo_duplex_fingerprint_data_WT_TCR.tsv | cut -f5 | sort -u > ../data/peptides.txt

/home/ctools/netMHCpan-4.2/netMHCpan -p ../data/peptides.txt -a HLA-B08:01 > ../results/netmhcpan_baseline.txt
