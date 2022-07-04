#!/usr/bin/env bash
set -euo pipefail
rg -t r '".*/.*"'  -o -I -N  | grep -v "{odir}\|.pdf\|.xlsx\|utils\|=\|temp/\|mkdir\|>\|module"  | sed 's/"//g' > files_read_by_r_data.txt 


