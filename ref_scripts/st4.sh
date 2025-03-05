#!/bin/bash

CHR=$1
CHUNK=$2

while IFS="" read -r LINE || [ -n "$LINE" ];
   do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)

        ./GLIMPSE2_split_reference_static \
            --reference /mnt/project/data/hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.bcf \
            --map chr${CHR}.b38.gmap.gz --input-region ${IRG} --output-region ${ORG} --output \
            HGDP_ref.chr${CHR}_ref  --threads 32
   done < ./$CHUNK


