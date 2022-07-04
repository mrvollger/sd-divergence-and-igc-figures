#!/usr/bin/env bash
set -eu 


zcat data/sd-divergence-results/all_snv_exploded.bed.gz \
  | grep -v "CHM1_2\|GRCh38_2" \
  | csvtk -tT -C '$' cut -f "#CHROM",POS,END,ID,anno_SD,anno_Unique \
  | bedtools sort \
  | uniq -c  \
  | awk -v OFS=$'\t' '{print $2,$3,$4,$5,$6,$7,$1}' \
  | bedtools intersect -a - \
    -b data/haplotype_coverage_bed_graph.bed \
    -wa -wb \
  | bgzip \
  > data/sd-divergence-results/snv_counts.bed.gz


zcat data/sd-divergence-results/snv_counts.bed.gz | head -n 10

echo "done"

exit
  #| head -n 1000000 \


