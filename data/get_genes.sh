#!/usr/bin/env bash
set -euo pipefail
wget https://t2t.gi.ucsc.edu/chm13/dev/t2t-chm13-v1.1/CAT_V4/CHM13.combined.v4.bb
bigBedToBed -header CHM13.combined.v4.bb CHM13.combined.v4.bb.bed 
bedtools bed12tobed6 -i CHM13.combined.v4.bb.bed  > CHM13.combined.v4.bb.bed6
