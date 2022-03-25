#!/usr/bin/env bash
set -euo pipefail



zcat data/gene-conversion/merged_acceptor.bed.gz \
  | grep -v "chrY\|CHM1_2\|GRCh38_2" \
  | bedtools sort -i - -g data/anno/chm13_v1.1_plus38Y.fasta.fai \
  | bedtools genomecov -bg -g data/anno/chm13_v1.1_plus38Y.fasta.fai -i - > data/gene-conversion/merged_acceptor.bed_graph.bed


