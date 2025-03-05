#!/bin/bash

#set this to the exome sequence directory that you want (should contain PLINK formatted files)
input_file_dir="/data/"
data_file_dir="/testrun/"
ref_file_dir="/refs/"


# Normaile making all multialleelic site biallelic
# extract site information

for CHR in {1..22}; do

run_cmd="wget  https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_chunk_static; \
          chmod a+x GLIMPSE2_chunk_static; \
          ./GLIMPSE2_chunk_static --input  /mnt/project/testrun/HGDP_ref.chr${CHR}.sites.vcf.gz --region chr${CHR} --sequential \
              --window-mb 6 --output chunks.chr${CHR}.txt --map chr${CHR}.b38.gmap.gz --threads 34 "

         dx run swiss-army-knife -iin="${ref_file_dir}maps/genetic_maps.b38/chr${CHR}.b38.gmap.gz" \
             -icmd="${run_cmd}" --tag="norm-extractsites" --instance-type "mem1_ssd2_v2_x36" \
             --destination="${project}:${data_file_dir}" --brief --yes

done
