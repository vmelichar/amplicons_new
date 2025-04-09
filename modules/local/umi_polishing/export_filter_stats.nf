process EXPORT_FILTER_STATS {
    conda '/home/melichv/miniconda3/envs/amplicons'
    publishDir "${params.output}/${sample}/export", pattern: "*.csv", mode: 'copy'

    input:
        path export_stats
        tuple val ( sample ), val ( hs_index ), val ( low_clusters_counts )
    
    output:
        tuple val( "${sample}" ), path ( "*.csv" )

    script:
        def counts = low_clusters_counts.join(' ')
        def hs_idx = hs_index.join(' ')

    """
        python ${export_stats} -i "${workflow.launchDir}/${params.output}" -l ${counts} -x ${hs_idx}
    """
}