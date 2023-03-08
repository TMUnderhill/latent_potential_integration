#!/bin/bash
samples=("fap_bat" "fap_brain" "fap_cerebellum" "fap_heart" "fap_kidney" "fap_lung" "fap_marrow" "fap_muscle" "fap_pancreas" "fap_salglnd" "fap_spleen" "fap_wat" "hsc_liver" "lsec_liver" "peri_bat" "peri_brain" "peri_cerebellum" "peri_heart" "peri_kidney" "peri_marrow" "peri_muscle" "peri_pancreas" "peri_salglnd" "peri_spleen" "peri_wat")

for s in ${samples[@]}; do
        tissue=$(echo $s | cut -d'_' -f 2)
        celltype=$(echo $s | cut -d'_' -f 1)
        echo $tissue ${celltype^^}
        nohup ./pipeline_final.R --rdata $tissue.Rdata --celltype ${celltype^^} --tssranges wstssranges.Rdata --fragments "fragments_"$tissue".tsv.gz" --overlapping overlapingTSS500.csv --nbatches 10 --cpm_cutoff 1 &
done