#!/usr/bin/env bash
set -euo pipefail
chr=${1:-chr20}; map=$2; out=${3:-ceph_$chr}; here=$(cd "$(dirname "$0")" && pwd); bin=${GIBD_BIN:-$here/../..}
mkdir -p "$out/tmp" && cd "$out"
bcftools view -r "$chr" -Oz -o truth.vcf.gz https://platinum-pedigree-data.s3.amazonaws.com/variants/small_variant_truthset/GRCh38/CEPH1463.GRCh38.family-truthset.ov.vcf.gz
bcftools index -t truth.vcf.gz
curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$chr.fa.gz" | gzip -dc > ref.fa
vg autoindex --workflow giraffe -r ref.fa -v truth.vcf.gz -p ped -T tmp
vg convert -f ped.giraffe.gbz > ped.gfa
awk -F'\t' '$1=="S"{print $2"\t"length($3)}' ped.gfa > ped.nlen
"$bin/gfa2gbwt" ped.gfa ped
"$bin/gibd" --nlen ped.nlen --ref-prefix "$chr" --map "$map" -o ped.ibd.tsv ped.gbwt
python3 "$here/score_pedigree.py" ped.ibd.tsv --ped "$here/CEPH1463.ped" --iv "$here/ceph.GRCh38.viterbi.oa.csv.gz" --map "$map" --chrom "$chr" --out score.tsv
