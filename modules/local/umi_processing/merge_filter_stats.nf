process MERGE_FILTER_STATS {
    tag "${sample}_${target}"
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.tsv", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( stats_file )
        val ( type )
        path merge_filter_stats_python
    
    output:
        path "*stats.tsv"

    script:
        def write_report = params.write_reports ? "--tsv" : ""

    """
        python ${merge_filter_stats_python} \
        --filter_tsv ${stats_file} \
        $write_report \
        -o .
    """
}