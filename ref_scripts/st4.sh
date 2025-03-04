#!/bin/bash

mkdir -p reference_panel/split

for i in {1..22}; do


while IFS="" read -r LINE || [ -n "$LINE" ];
   do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)

        GLIMPSE2_split_reference \ 
            --reference ~/ariel_research/lp-wgs_research/hgdp_ref/phased_haplotypes_v2/hgdp1kgp_chr${i}.filtered.SNV_INDEL.phased.shapeit5.bcf \
            --map ~/GLIMPSE/maps/genetic_maps.b38/chr${i}.b38.gmap.gz --input-region ${IRG} --output-region ${ORG} --output \
            reference_panel/split/HGDP_ref.chr${i}_ref  --threads 32
   done < reference_panel/chunks.chr${i}.txt



done
