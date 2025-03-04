#!/bin/bash


for CHR in {1..22}; do

    ~/GLIMPSE/bin/GLIMPSE2_chunk --input  ./reference_panel/HGDP_ref.chr${CHR}.sites.vcf.gz --region chr${CHR} --sequential --window-mb 6 --output \
         ./reference_panel/chunks.chr${CHR}.txt --map ~/GLIMPSE/maps/genetic_maps.b38/chr${CHR}.b38.gmap.gz --threads 32

done
