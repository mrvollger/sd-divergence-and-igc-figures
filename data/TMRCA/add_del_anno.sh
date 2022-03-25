#!/usr/bin/env bash
set -euo pipefail


ls ./sd.age.txt ./uniq.age.txt \
  | parallel -n 1 'bedtools annotate -i {} -files large_clint_dels.bed | bedtools nuc -fi ../anno/chm13_v1.1_plus38Y.fasta -bed - | grep -v "^#" > del_anno_{/};'


exit

ls 10*bed \
  | parallel -n 1 'bedtools annotate -i {} -files large_clint_dels.bed > del_anno_{/};'

awk '$5>0.2' ./del_anno_10_kbp_windows_with_50_coverage_and_90_percent_sd.bed > bad_sd_windows.bed
awk '$4>0.2' ./del_anno_10_kbp_windows_with_50_coverage_and_unique.10kWindows.bed > bad_unique_windows.bed
