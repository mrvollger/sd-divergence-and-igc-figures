#!/usr/bin/env bash
set -euo pipefail


bedtools annotate \
  -counts \
  -i regions_to_validate.bed \
  -files ../gene-conversion/tables/*acceptor.bed \
  -names $(ls ../gene-conversion/tables/*acceptor.bed | sed 's#../gene-conversion/tables/##g' | sed 's/.acceptor.bed//g') \
  > regions_to_validate.with.igc.counts.bed

