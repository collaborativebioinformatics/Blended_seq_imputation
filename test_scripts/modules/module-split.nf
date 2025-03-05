// Module for splitting reference into binary chunks
process SPLIT_REFERENCE {
    tag "split_${chr}_chunk${chunk_id}"
    container params.docker_image
    publishDir "${params.output_dir}/reference_panels", mode: params.publish_mode
    
    input:
    tuple val(chr), val(chunk_id), val(region), path(bcf_file), path(bcf_index)
    
    output:
    tuple val(chr), val(chunk_id), path("${params.sample_name}.${chr}.${chunk_id}.bin")
    
    script:
    """
    GLIMPSE2_split_reference --input ${bcf_file} \
        --region ${region} \
        --output ${params.sample_name}.${chr}.${chunk_id}.bin
    """
}

// Function to apply the process
def splitReference(input_channel) {
    return SPLIT_REFERENCE(input_channel)
}
