// Module for extracting site information
process EXTRACT_SITE_INFO {
    tag "extract_${chr}"
    container params.docker_image
    publishDir "${params.output_dir}/sites", mode: params.publish_mode
    
    input:
    tuple val(chr), path(bcf_file), path(bcf_index)
    
    output:
    tuple val(chr), path("${chr}.sites.vcf.gz"), path("${chr}.sites.vcf.gz.tbi")
    
    script:
    """
    bcftools view -G -m2 -M2 -v snps,indels -Oz -o ${chr}.sites.vcf.gz ${bcf_file}
    tabix -p vcf ${chr}.sites.vcf.gz
    """
}

// Function to apply the process
def extractSiteInfo(input_channel) {
    return EXTRACT_SITE_INFO(input_channel)
}
