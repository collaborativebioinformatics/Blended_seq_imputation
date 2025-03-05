#!/usr/bin/env nextflow

// Modular Imputation Reference Panel Pipeline
// Date: 2025-03-04

// Enable DSL2 for modularity
nextflow.enable.dsl=2

// Default parameters
params {
    // Input parameters with sensible defaults
    input_pattern = null               // Required: Input pattern for BCF files (e.g., "/data/samples_chr*.bcf")
    input_index_pattern = null         // Optional: Input pattern for index files (default: auto-detect)
    sample_name = "reference_panel"    // Default sample name for output files
    
    // Chromosome handling
    chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" // Default chromosomes to process
    chrX_mode = "standard"             // Options: standard, split (for handling PAR/non-PAR regions)
    chrX_non_par_pattern = null        // For split mode: pattern for X non-PAR
    chrX_par1_pattern = null           // For split mode: pattern for X PAR1
    chrX_par2_pattern = null           // For split mode: pattern for X PAR2
    
    // Processing parameters
    skip_multiallelic = false          // Skip converting multiallelic sites if input is already biallelic
    sites_only = false                 // Skip site extraction if sites VCF already exists
    
    // Output parameters
    output_dir = "results"             // Output directory
    publish_mode = "copy"              // Mode for publishDir (copy, link, symlink, move)
    
    // GLIMPSE2 parameters
    window_size = 2000000              // Window size for chunking (2Mb default)
    buffer_size = 200000               // Buffer size for chunks (200kb default)
    
    // Execution parameters
    docker_image = "pjgreer/glimpse2plus:v0.1"  // Docker image to use
    max_cpus = 16                      // Maximum CPUs per process
    max_memory = "32 GB"               // Maximum memory per process
    max_time = "24h"                   // Maximum time per process
}

// Validate required parameters
if (params.input_pattern == null) {
    error "Input pattern must be specified using --input_pattern parameter"
}

// Process input parameters 
// Parse chromosome list
def chrList = params.chromosomes.split(',')

// Print workflow header with parameters
def printHeader() {
    log.info """
    ==============================================
    MODULAR BGE IMPUTATION REFERENCE PANEL PIPELINE
    ==============================================
    input_pattern        : ${params.input_pattern}
    input_index_pattern  : ${params.input_index_pattern ?: "auto-detected"}
    sample_name          : ${params.sample_name}
    chromosomes          : ${params.chromosomes}
    chrX_mode            : ${params.chrX_mode}
    skip_multiallelic    : ${params.skip_multiallelic}
    sites_only           : ${params.sites_only}
    output_dir           : ${params.output_dir}
    window_size          : ${params.window_size}
    buffer_size          : ${params.buffer_size}
    docker_image         : ${params.docker_image}
    """
}

// Print header
printHeader()

// Define modules

// Module: Convert multiallelic to biallelic
process convertMultiallelicToBiallelic {
    tag "chr${chr}"
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

// Module: Extract site information
process extractSiteInfo {
    tag "chr${chr}"
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

// Module: Create chunks using GLIMPSE2_chunk
process createChunks {
    tag "chr${chr}"
    container params.docker_image
    publishDir "${params.output_dir}/chunks", mode: params.publish_mode
    
    input:
    tuple val(chr), path(sites_vcf), path(sites_vcf_index)
    
    output:
    tuple val(chr), path("${chr}.chunks.txt")
    
    script:
    def chr_id = chr == "X" ? "X" : chr
    """
    GLIMPSE2_chunk --input ${sites_vcf} --region ${chr_id} --window-size ${params.window_size} --buffer-size ${params.buffer_size} --output ${chr}.chunks.txt
    """
}

// Module: Split reference into binary chunks
process splitReferenceIntoChunks {
    tag "chr${chr}_chunk${chunk_id}"
    container params.docker_image
    publishDir "${params.output_dir}/reference_panels", mode: params.publish_mode
    
    input:
    tuple val(chr), val(chunk_id), val(region), path(bcf_file), path(bcf_index)
    
    output:
    tuple val(chr), val(chunk_id), path("${params.sample_name}.${chr}.${chunk_id}.bin")
    
    script:
    """
    GLIMPSE2_split_reference --input ${bcf_file} --region ${region} --output ${params.sample_name}.${chr}.${chunk_id}.bin
    """
}

// Module: Summarize results
process summarizeResults {
    tag "summary"
    container params.docker_image
    publishDir "${params.output_dir}", mode: params.publish_mode
    
    input:
    path('panels/*')
    
    output:
    path("reference_panel_summary.txt")
    
    script:
    """
    echo "Reference Panel Summary" > reference_panel_summary.txt
    echo "======================" >> reference_panel_summary.txt
    echo "Total panels: \$(ls panels/ | wc -l)" >> reference_panel_summary.txt
    echo "Panels by chromosome:" >> reference_panel_summary.txt
    for chr in \$(ls panels/ | cut -d'.' -f2 | sort -u); do
        count=\$(ls panels/ | grep ".\${chr}." | wc -l)
        echo "  Chromosome \${chr}: \${count} panels" >> reference_panel_summary.txt
    done
    echo "======================" >> reference_panel_summary.txt
    echo "Panel details:" >> reference_panel_summary.txt
    ls -lh panels/ >> reference_panel_summary.txt
    """
}

// Sub-workflow: Process single chromosome
workflow processChr {
    take:
        chr_channel
    
    main:
        // Apply appropriate workflow based on parameters
        if (params.skip_multiallelic) {
            biallelic_files = chr_channel
        } else {
            biallelic_files = convertMultiallelicToBiallelic(chr_channel)
        }
        
        if (params.sites_only) {
            site_files = chr_channel.map { chr, bcf, index -> 
                // Assume sites files are named as ${chr}.sites.vcf.gz
                return tuple(chr, file("${params.output_dir}/sites/${chr}.sites.vcf.gz"), file("${params.output_dir}/sites/${chr}.sites.vcf.gz.tbi"))
            }
        } else {
            site_files = extractSiteInfo(biallelic_files)
        }
        
        chunk_files = createChunks(site_files)
        
        // Parse chunks file to get regions
        chunk_regions = chunk_files.flatMap { chr, chunks_file ->
            def lines = chunks_file.readLines()
            def result = []
            
            for (int i = 1; i < lines.size(); i++) {  // Skip header line
                def fields = lines[i].split()
                // chunk_id, contig, physical_pos_start, physical_pos_end, ...
                def chunk_id = fields[0]
                def region = fields[2]
                result << tuple(chr, chunk_id, region)
            }
            
            return result
        }
        
        // Combine with biallelic files to prepare for splitting
        split_inputs = chunk_regions.combine(biallelic_files, by: 0)
            .map { chr, chunk_id, region, bcf, index -> 
                tuple(chr, chunk_id, region, bcf, index)
            }
        
        // Split reference
        reference_panels = splitReferenceIntoChunks(split_inputs)
    
    emit:
        reference_panels
}

// Main workflow
workflow {
    // Create channel for input files
    if (params.input_index_pattern) {
        // If index pattern is provided
        Channel.fromFilePairs("${params.input_pattern};${params.input_index_pattern}", size: 2)
            .map { pattern, files -> 
                def chr = pattern.replaceAll(/.*chr/, "").replaceAll(/\..*/, "")
                return tuple(chr, files[0], files[1])
            }
            .filter { chr, bcf, index -> chr in chrList }
            .set { input_files }
    } else {
        // Auto-detect index files
        Channel.fromPath(params.input_pattern)
            .map { bcf -> 
                def chr = bcf.name.toString().replaceAll(/.*chr/, "").replaceAll(/\..*/, "")
                def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : 
                            file("${bcf}.tbi").exists() ? file("${bcf}.tbi") : null
                if (index == null) {
                    error "Cannot find index for ${bcf}. Please provide index or generate it."
                }
                return tuple(chr, bcf, index)
            }
            .filter { chr, bcf, index -> chr in chrList }
            .set { input_files }
    }
    
    // Handle X chromosome in split mode
    if (params.chrX_mode == "split" && chrList.contains("X")) {
        // Remove X from standard processing
        input_files = input_files.filter { chr, bcf, index -> chr != "X" }
        
        // Create channels for X chromosome parts
        if (params.chrX_non_par_pattern && params.chrX_par1_pattern && params.chrX_par2_pattern) {
            // Create channel for X non-PAR
            Channel.fromPath(params.chrX_non_par_pattern)
                .map { bcf -> 
                    def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : 
                                file("${bcf}.tbi").exists() ? file("${bcf}.tbi") : null
                    return tuple("X_non_par", bcf, index)
                }
                .set { chrX_non_par }
                
            // Create channel for X PAR1
            Channel.fromPath(params.chrX_par1_pattern)
                .map { bcf -> 
                    def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : 
                                file("${bcf}.tbi").exists() ? file("${bcf}.tbi") : null
                    return tuple("X_par1", bcf, index)
                }
                .set { chrX_par1 }
                
            // Create channel for X PAR2
            Channel.fromPath(params.chrX_par2_pattern)
                .map { bcf -> 
                    def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : 
                                file("${bcf}.tbi").exists() ? file("${bcf}.tbi") : null
                    return tuple("X_par2", bcf, index)
                }
                .set { chrX_par2 }
                
            // Combine X parts with autosomal inputs
            input_files
                .mix(chrX_non_par)
                .mix(chrX_par1)
                .mix(chrX_par2)
                .set { all_inputs }
        } else {
            log.warn "X chromosome split mode enabled but patterns not provided, using standard mode for X"
            all_inputs = input_files
        }
    } else {
        all_inputs = input_files
    }
    
    // Process each chromosome
    reference_panels = processChr(all_inputs)
    
    // Collect all panels for summary
    reference_panels.map { chr, chunk_id, panel -> panel }
        .collect()
        .set { all_panels }
    
    // Generate summary
    summarizeResults(all_panels.collect())
}