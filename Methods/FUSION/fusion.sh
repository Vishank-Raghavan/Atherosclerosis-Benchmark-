#!/usr/bin/env bash
# Set the path to the FUSION R script
FUSION_SCRIPT="./fusion_twas-master/FUSION.assoc_test.R"

for CHR in {1..22}; do
  echo "running $CHR"
  Rscript $FUSION_SCRIPT \
    --sumstats ./CAD.sumstats \
    --weights ./GTExv8.ALL.Artery_Coronary.pos \
    --weights_dir ./ \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr $CHR \
    --out FUSION_TWAS_Artery_Coronary_results_chr${CHR}.dat \
    --GWASN 184305
done
