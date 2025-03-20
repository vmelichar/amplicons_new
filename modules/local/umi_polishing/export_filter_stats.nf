process EXPORT_FILTER_STATS {
    conda '/home/melichv/miniconda3/envs/amplicons'
    publishDir "${params.output}/${sample}/export", pattern: "*.csv", mode: 'copy'

    input:
        path export_stats
        val ( dummy )
    
    output:
        tuple val( "${sample}" ), path ( "*.csv" )

    """
        python ${export_stats} -i "${workflow.launchDir}/${params.output}"
    """
}