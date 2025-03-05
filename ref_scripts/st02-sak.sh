#!/bin/bash

#set this to the exome sequence directory that you want (should contain PLINK formatted files)
input_file_dir="/data/"
data_file_dir="/testrun/"
ref_file_dir="/refs/"


# Normaile making all multialleelic site biallelic
# extract site information

for CHR in {1..22}; do

run_cmd="bcftools norm -m -any /mnt/project/data/hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.bcf \
             -Ob --threads 4 > hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.norm.bcf; \
         bcftools index -f hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.norm.bcf; \
         bcftools view -G -Oz -o HGDP_ref.chr${CHR}.sites.vcf.gz hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.norm.bcf; \
         bcftools index -f HGDP_ref.chr${CHR}.sites.vcf.gz "


         dx run swiss-army-knife -iin="${ref_file_dir}maps/genetic_maps.b38/chr${CHR}.b38.gmap.gz" \
             -icmd="${run_cmd}" --tag="norm-extractsites" --instance-type "mem1_ssd2_v2_x8" \
             --destination="${project}:${data_file_dir}" --brief --yes

done
