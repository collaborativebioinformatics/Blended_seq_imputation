// Module for creating a summary of the reference panels
process CREATE_SUMMARY {
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

// Function to apply the process
def createSummary(input_channel) {
    return CREATE_SUMMARY(input_channel)
}
