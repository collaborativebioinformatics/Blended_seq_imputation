// Module for creating chunks using GLIMPSE2_chunk
process CREATE_CHUNKS {
    tag "chunk_${chr}"
    container params.docker_image
    publishDir "${params.output_dir}/chunks", mode: params.publish_mode
    
    input:
    tuple val(chr), path(sites_vcf), path(sites_vcf_index)
    
    output:
    tuple val(chr), path("${chr}.chunks.txt")
    
    script:
    def chr_id = chr.startsWith("X") ? "X" : chr
    """
    GLIMPSE2_chunk --input ${sites_vcf} --region ${chr_id} \
        --window-size ${params.window_size} \
        --buffer-size ${params.buffer_size} \
        --output ${chr}.chunks.txt
    """
}

// Function to apply the process
def createChunks(input_channel) {
    return CREATE_CHUNKS(input_channel)
}
