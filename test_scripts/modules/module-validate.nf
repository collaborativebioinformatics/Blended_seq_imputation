// Module to validate and prepare input files
process VALIDATE_INPUTS {
    tag "validate_${chr}"
    container params.docker_image
    
    input:
    tuple val(chr), path(bcf_file), path(bcf_index)
    
    output:
    tuple val(chr), path("${chr}.validated.bcf"), path("${chr}.validated.bcf.csi")
    
    script:
    """
    # Simple validation: check if file is readable
    bcftools view -h ${bcf_file} > /dev/null
    
    # Create a symlink with standard name
    ln -s ${bcf_file} ${chr}.validated.bcf
    ln -s ${bcf_index} ${chr}.validated.bcf.csi
    """
}

// Function to create input channel from file patterns
def validateInputs(input_pattern, index_pattern, chrList) {
    if (input_pattern == null) {
        error "Input pattern must be specified using --input_pattern parameter"
    }
    
    if (index_pattern) {
        // If index pattern is provided
        return Channel.fromFilePairs("${input_pattern};${index_pattern}", size: 2)
            .map { pattern, files -> 
                def chr = pattern.replaceAll(/.*chr/, "").replaceAll(/\..*/, "")
                return tuple(chr, files[0], files[1])
            }
            .filter { chr, bcf, index -> chr in chrList }
    } else {
        // Auto-detect index files
        return Channel.fromPath(input_pattern)
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
    }
}
