#!/usr/bin/env bash
set -euo pipefail
IGC="../data/gene-conversion/acceptor_gene_conversion_windows.bed.gz"
FAI="../../../chm13_t2t/assemblies/chm13_v1.1_plus38Y.fasta.fai"
COV="../data/haplotype_coverage_bed_graph.bed"

W=1000
zcat $IGC | tail -n +2 \
  | grep -wv "GRCh38_2\|CHM1_2\|chrY" \
  | bedtools sort -i - \
  | cut -f 11-13 > d.bed
zcat $IGC | tail -n +2 \
  | grep -wv "GRCh38_2\|CHM1_2\|chrY" \
  | bedtools sort -i - \
  | cut -f 1-3 > a.bed

rb bl -r a.bed 
rb bl -r d.bed

bedtools intersect -f 0.5 -C \
  -a <(bedtools makewindows -w $W \
    -b <(cat a.bed d.bed | bedtools sort -i - | bedtools slop -i - -b $W -g $FAI | bedtools merge -i - ) \
  ) \
  -b a.bed \
  > a.2.bed

bedtools intersect -f 0.5 -C \
  -a <(bedtools makewindows -w $W \
    -b <(cat a.bed d.bed | bedtools sort -i - | bedtools slop -i - -b $W -g $FAI | bedtools merge -i - ) \
  ) \
  -b d.bed \
  > d.2.bed

printf "#chr\tstart\tend\tacceptor\tdonor\tcov_chr\tcov_start\tcov_end\thap_cov\n" \
  > ../data/windowed_acceptor_donor_counts.bed

paste a.2.bed <(cut -f 4 d.2.bed) \
  | bedtools intersect -a - -b $COV -wa -wb \
  >> ../data/windowed_acceptor_donor_counts.bed

exit
bedtools coverage -counts \
