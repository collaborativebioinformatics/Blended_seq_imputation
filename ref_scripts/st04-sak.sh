#!/bin/bash

#set this to the exome sequence directory that you want (should contain PLINK formatted files)
input_file_dir="/data/"
data_file_dir="/testrun/"
ref_file_dir="/refs/"


# Normaile making all multialleelic site biallelic
# extract site information

for CHR in {1..19}; do

run_cmd="wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_split_reference_static; \
           chmod a+x GLIMPSE2_split_reference_static; \
           cp /mnt/project/scripts/st4.sh . ; \
           bash st4.sh ${CHR} chunks.chr${CHR}.txt"

         dx run swiss-army-knife -iin="${ref_file_dir}maps/genetic_maps.b38/chr${CHR}.b38.gmap.gz" \
             -iin="${data_file_dir}/chunks.chr${CHR}.txt" \
             -icmd="${run_cmd}" --tag="norm-extractsites" --instance-type "mem1_ssd2_v2_x36" \
             --destination="${project}:${data_file_dir}" --brief --yes

done
