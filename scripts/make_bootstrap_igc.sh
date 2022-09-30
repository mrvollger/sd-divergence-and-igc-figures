#!/usr/bin/env bash
set -euo pipefail

(printf "#chrm\tst\ten\tbootstrap\n"; seq $1 \
  | parallel -n 1 \
  'bedtools shuffle -maxTries 1000000 -i <(zcat data/gene-conversion/merged_acceptor.bed.gz  | cut -f 1-3) -incl <(bedtools merge -i data/anno/chm13_v1.1_plus38Y.SDs.bed) -g data/anno/chm13_v1.1_plus38Y.fasta.fai | sed "s/$/\t{}/g" ' \
  | bedtools sort -i - ) \
  | bgzip -@ 30 





