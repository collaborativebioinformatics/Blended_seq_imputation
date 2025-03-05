// Module for converting multiallelic sites to biallelic sites
process CONVERT_MULTIALLELIC {
    tag "convert_${chr}"
    container params.docker_image
    publishDir "${params.output_dir}/biallelic", mode: params.publish_mode
    
    input:
    tuple val(chr), path(bcf_file), path(bcf_index)
    
    output:
    tuple val(chr), path("${chr}.biallelic.bcf"), path("${chr}.biallelic.bcf.csi")
    
    script:
    """
    bcftools norm -m-any -Ob -o ${chr}.biallelic.bcf ${bcf_file}
    bcftools index ${chr}.biallelic.bcf
    """
}

// Function to apply the process
def convertMultiallelic(input_channel) {
    return CONVERT_MULTIALLELIC(input_channel)
}
