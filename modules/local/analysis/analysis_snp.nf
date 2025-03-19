process ANALYSIS_SNP {
    // publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.tsv", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( det_umi ), path ( extr_syn ), path ( extr_umi )
        path positions_file
        path analysis_script
    
    output:
        // path "*stats.tsv", emit: stats_tsv

    script:

    """
        python ${analysis_script} \
        
    """
}