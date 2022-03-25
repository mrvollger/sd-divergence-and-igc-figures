#!/usr/bin/env bash
set -euo pipefail

mkdir -p temp
rm -rf temp/*.bed

if [ "x" == "y" ]; then
  S3="s3://human-pangenomics/submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/JAN_09_2022/FINAL_HIFI_BASED/FLAGGER_HIFI_PROJECTED_SIMPLIFIED_BEDS/ALL/"
  S3="s3://human-pangenomics/submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/JAN_09_2022/FINAL_HIFI_BASED/FLAGGER_HIFI_ASM_SIMPLIFIED_BEDS/ALL/"
  aws s3 --no-sign-request sync $S3 temp/.
  cat temp/*.bed | grep -v "^track" | awk '$4!="Hap"' | bgzip > all.issues.bed.gz
  
  zcat all.issues.bed.gz | cut -f 1 | sed 's/#/_/g' | sed 's/_J.*//g' | sort | uniq  > samples.txt
fi

if [ "x" == "y" ]; then
  bedtools merge -i ../anno/chm13_v1.1_plus38Y.SDs.bed > tmp.sds.bed 
  ls ~/assembly_breaks/asm-to-reference-alignment/results/CHM13_V1.1_v2.24_asm20/paf_trim_and_break/*paf \
    | parallel -n 1 \
      $'awk \' $4-$3 > 1000000 \' {} | rb -v liftover --bed tmp.sds.bed ' \
    > sds.paf
fi 

cat samples.txt \
  | parallel -n 1 \
  'zcat ../callable/{}*bed.gz | csvtk cut -tT -f 4,5,6,1,2,3 ' \
  | bgzip -@ 30 > callable.bed.gz


bedtools intersect -wb -a callable.bed.gz -b all.issues.bed.gz > issues.in.callable.bed 

bedtools intersect -wb -a <(cut -f 1,3,4 sds.paf ) -b issues.in.callable.bed > issues.in.callable.sds.bed

rb bl -r issues.in.callable.bed
rb bl -r issues.in.callable.sds.bed
rb bl -c 4 -r <( cut -f 1,2,3,13 issues.in.callable.sds.bed)


rm -rf temp/*.bed
rm tmp.*.bed

