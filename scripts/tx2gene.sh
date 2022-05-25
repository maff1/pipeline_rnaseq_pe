#!/bin/bash
cd /mnt/marcodata/PROJECTS/snakemake_uzh/reference
gtf=Homo_sapiens.GRCh38.105.gtf
test=$(zless -S $gtf | grep -v "#" | awk '$3=="transcript"' | head -n 1| cut -f9 | tr -s ";" " " | awk '{print$3}' | sort | uniq | sed 's/"//g')
if [[ $test == "transcript_id" ]]; then
zless -S $gtf | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort | uniq | sed 's/"//g' > txp2gene.tsv
elif [[ $test == "gene_version" ]]; then
echo "Separate version field (ensembl, non-gencode transcriptome, eg. rat, etc)"
zless -S $gtf | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$6 "." $8"\t"$2 "." $4}' | sort | uniq | sed 's/"//g' > txp2gene.tsv
fi