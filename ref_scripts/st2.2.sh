# Multiple processing commands piped together

for CHR in {20,21,22,}; do

# Normaile making all multialleelic site biallelic
 bcftools norm -m -any ~/ariel_research/lp-wgs_research/hgdp_ref/phased_haplotypes_v2/hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.bcf \
     -Ob --threads 4 > hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.norm.bcf 
 bcftools index -f hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.norm.bcf

# extract site information
 bcftools view -G -Oz -o reference_panel/HGDP_ref.chr${CHR}.sites.vcf.gz hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.norm.bcf
 bcftools index -f reference_panel/HGDP_ref.chr${CHR}.sites.vcf.gz

done
