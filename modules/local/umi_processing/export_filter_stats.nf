process EXPORT_FILTER_STATS {
    publishDir "${params.output}/${sample}/export", pattern: "*.csv", mode: 'copy'

    input:
        path export_stats
        path cluster_stats
    
    output:
        tuple val( "${sample}" ), path ( "*.csv" )

    """
        python ${export_stats} -i ${params.output}
    """
}