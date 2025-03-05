#!/usr/bin/env nextflow

// Pipeline parameters with default values and flexibility
params.input_bcf_list = "samples.txt"
params.output_dir = "reference_panel"
params.map_dir = "maps"
params.threads = 4

process filter_bcf {
    publishDir "${params.output_dir}", mode: "copy"

    input:
        path input_bcf

    output:
        tuple path("${input_bcf.getSimpleName()}.qc.bcf"), path("${input_bcf.getSimpleName()}.qc.bcf.csi")

    script:
    """
    bcftools norm -m -any ${input_bcf} -Ou --threads ${params.threads} | 
    bcftools view -m 2 -M 2 -v snps,indels --threads ${params.threads} -Ob -o ${input_bcf.getSimpleName()}.qc.bcf
    bcftools index -f ${input_bcf.getSimpleName()}.qc.bcf --threads ${params.threads}
    """
}

process extract_sites {
    publishDir "${params.output_dir}", mode: "copy"

    input:
        tuple path(filtered_bcf), path(filtered_bcf_index)

    output:
        tuple path("${filtered_bcf.baseName}.sites.vcf.gz"), path("${filtered_bcf.baseName}.sites.vcf.gz.tbi")

    script:
    """
    bcftools view -G -Oz -o ${filtered_bcf.baseName}.sites.vcf.gz ${filtered_bcf}
    tabix -p vcf ${filtered_bcf.baseName}.sites.vcf.gz
    """
}

process chunks {
    publishDir "${params.output_dir}/chunks", mode: "copy"

    input:
        tuple path(sites_vcf), path(sites_vcf_index)
        path(map)

    output:
        path "chunks.${sites_vcf.getSimpleName()}.txt"

    script:
    def filename = sites_vcf.getSimpleName()
    def matcher = (filename =~ /(chr\d+)/)
    def chr = matcher[0][1]

    """
    GLIMPSE2_chunk --input ${sites_vcf} --region ${chr} --sequential --window-mb 6 \
    --output chunks.${sites_vcf.getSimpleName()}.txt \
    --map ${chr}.b38.gmap.gz --threads ${params.threads}
    """
}

process split_reference {
    publishDir "${params.output_dir}/split", mode: "copy"

    input:
        tuple path(filtered_bcf), path(filtered_bcf_index)
        path(chunk_file)
        path(map)

    output:
        path "*"

    script:
    def filename = filtered_bcf.getSimpleName()
    def matcher = (filename =~ /(chr\d+)/)
    def chr = matcher[0][1]

    """
    while IFS="" read -r LINE || [ -n "\$LINE" ];
    do
        printf -v ID "%02d" \$(echo \$LINE | cut -d" " -f1)
        IRG=\$(echo \$LINE | cut -d" " -f3)
        ORG=\$(echo \$LINE | cut -d" " -f4)

        GLIMPSE2_split_reference \\
            --reference ${filtered_bcf} \\
            --map ${chr}.b38.gmap.gz \\
            --input-region \$IRG \\
            --output-region \$ORG \\
            --output HGDP_ref.${filtered_bcf.getSimpleName()}_ref \\
            --threads ${params.threads}
    done < chunks.${filtered_bcf.getSimpleName()}.txt
    """
}

workflow {
    // Create input channel for BCF files
    input_ch = Channel.fromPath(params.input_bcf_list)
                    .splitText()
                    .map { it.trim() }
                    .map { file(it) }
    
    map_ch = Channel.fromPath("${params.map_dir}/*.gmap.gz")
                    .collect()
                    .map { it }

    // Run the filter_bcf process for each input file
    filter_bcf_ch = filter_bcf(input_ch)

    // Run the extract_sites process for each filtered BCF file
    extract_sites_ch = extract_sites(filter_bcf_ch)

    // Run the chunks process for each sites VCF file
    chunks_ch = chunks(extract_sites_ch, map_ch).collect().map { it }

    // Run the split_reference process for each chunk file
    split_reference(filter_bcf_ch, chunks_ch, map_ch)
}
