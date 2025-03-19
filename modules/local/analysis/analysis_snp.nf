process ANALYSIS_SNP {
    publishDir "${params.output}/${sample}/${target}/analysis", pattern: "XXXX", mode: 'copy'

    input:
        tuple val( target ), val( sample ), path ( bam ), path ( bai ), path ( bed )
        path positions_file
        path variants_vcf
        path variants_tbi
        path analysis_script
    
    output:
        // path "*stats.tsv", emit: stats_tsv

    script:

    """
        python ${analysis_script} \
        ${bam} \
        ${bai} \
        ${variants_vcf}
        ${varinats_tbi}
        ${positions_file}
        .
        
    """
}