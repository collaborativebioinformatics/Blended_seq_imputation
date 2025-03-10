#!/usr/bin/env nextflow

// Modular BGE Imputation Reference Panel Pipeline
// Author: PJ Greer (modified)
// Date: 2025-03-04

// Enable DSL2 for modularity
nextflow.enable.dsl=2

// Import modules
include { convertMultiallelic } from './modules/convert_multiallelic'
include { extractSiteInfo } from './modules/extract_site_info'
include { createChunks } from './modules/create_chunks'
include { splitReference } from './modules/split_reference'
include { createSummary } from './modules/create_summary'

// Default parameters
params {
    // Input parameters with sensible defaults
    input_pattern = null               // Required: Input pattern for BCF files (e.g., "/data/samples_chr*.bcf")
    sample_name = "reference_panel"    // Default sample name for output files
    
    // Chromosome handling
    chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" // Default chromosomes to process
    chrX_mode = "standard"             // Options: standard, split (for handling PAR/non-PAR regions)
    chrX_non_par_pattern = null        // For split mode: pattern for X non-PAR
    chrX_par1_pattern = null           // For split mode: pattern for X PAR1
    chrX_par2_pattern = null           // For split mode: pattern for X PAR2
    
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
    threads = 4                        // Number of threads to use for each process
}

// Print workflow header with parameters
def printHeader() {
    log.info """
    ==============================================
    MODULAR BGE IMPUTATION REFERENCE PANEL PIPELINE
    ==============================================
    input_pattern        : ${params.input_pattern}
    sample_name          : ${params.sample_name}
    chromosomes          : ${params.chromosomes}
    chrX_mode            : ${params.chrX_mode}
    output_dir           : ${params.output_dir}
    window_size          : ${params.window_size}
    buffer_size          : ${params.buffer_size}
    docker_image         : ${params.docker_image}
    threads              : ${params.threads}
    """
}

// Main workflow
workflow {
    // Print header
    printHeader()
    
    // Parse chromosome list
    chrList = params.chromosomes.split(',')
    
    // Create input channel
    input_channel = createInputChannel(params.input_pattern, chrList)
    
    // Handle X chromosome in split mode if needed
    if (params.chrX_mode == "split" && chrList.contains("X")) {
        chrX_inputs = handleChrX(input_channel, chrList)
        input_channel = chrX_inputs
    }
    
    // Process through conversion
    biallelic_files = convertMultiallelic(input_channel)
    
    // Extract site information
    site_files = extractSiteInfo(biallelic_files)
    
    // Create chunks
    chunk_files = createChunks(site_files)
    
    // Parse chunk files to get regions
    chunk_regions = parseChunkFiles(chunk_files)
    
    // Combine with biallelic files to prepare for splitting
    split_inputs = chunk_regions.combine(biallelic_files, by: 0)
        .map { chr, chunk_id, region, bcf, index -> 
            tuple(chr, chunk_id, region, bcf, index)
        }
    
    // Split reference
    reference_panels = splitReference(split_inputs)
    
    // Collect all panels for summary
    all_panels = reference_panels.map { chr, chunk_id, panel -> panel }
        .collect()
    
    // Generate summary
    createSummary(all_panels.collect())
}

// Function to create input channel from file pattern
def createInputChannel(input_pattern, chrList) {
    if (input_pattern == null) {
        error "Input pattern must be specified using --input_pattern parameter"
    }
    
    return Channel.fromPath(input_pattern)
        .map { bcf -> 
            def chr = bcf.name.toString().replaceAll(/.*chr/, "").replaceAll(/\..*/, "")
            
            // Check if index exists, if not assume it has the same name with .csi extension
            def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : null
            
            // If index not found, try creating it using the error message
            if (index == null) {
                log.warn "Index for ${bcf} not found, will be created during processing"
                index = file("${bcf}.csi")
            }
            
            return tuple(chr, bcf, index)
        }
        .filter { chr, bcf, index -> chr in chrList }
}

// Function to handle X chromosome in split mode
def handleChrX(input_channel, chrList) {
    // Remove X from standard processing
    standard_channel = input_channel.filter { chr, bcf, index -> chr != "X" }
    
    // Create channels for X chromosome parts if patterns are provided
    if (params.chrX_non_par_pattern && params.chrX_par1_pattern && params.chrX_par2_pattern) {
        // Create channel for X non-PAR
        chrX_non_par = Channel.fromPath(params.chrX_non_par_pattern)
            .map { bcf -> 
                def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : null
                return tuple("X_non_par", bcf, index)
            }
            
        // Create channel for X PAR1
        chrX_par1 = Channel.fromPath(params.chrX_par1_pattern)
            .map { bcf -> 
                def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : null
                return tuple("X_par1", bcf, index)
            }
            
        // Create channel for X PAR2
        chrX_par2 = Channel.fromPath(params.chrX_par2_pattern)
            .map { bcf -> 
                def index = file("${bcf}.csi").exists() ? file("${bcf}.csi") : null
                return tuple("X_par2", bcf, index)
            }
            
        // Combine X parts with standard inputs
        return standard_channel
            .mix(chrX_non_par)
            .mix(chrX_par1)
            .mix(chrX_par2)
    } else {
        log.warn "X chromosome split mode enabled but patterns not provided, using standard mode for X"
        return input_channel
    }
}

// Function to parse chunk files
def parseChunkFiles(chunk_files) {
    return chunk_files.flatMap { chr, chunks_file ->
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
}