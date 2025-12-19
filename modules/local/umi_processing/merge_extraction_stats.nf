process MERGE_EXTRACTION_STATS {
    tag "${sample}_${target}"
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.tsv", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( det_umi ), path ( extr_syn ), path ( extr_umi )
        val ( type )
        path merge_extr_stats_python
    
    output:
        path "*stats.tsv", emit: stats_tsv
        tuple val( "${sample}" ), val ( "${target}" ), path( "extraction_synthetic_stats.tsv" ), emit: extr_synthetic_stats
        tuple val( "${sample}" ), val ( "${target}" ), path( "extraction_umi_stats.tsv" ), emit: extr_umi_stats
        val ( "${sample}" ), emit: dummy

    script:
        def write_report = params.write_reports ? "--tsv" : ""
        def cons = "${type}" == "consensus" ? "--cons" : ""

    """
        python ${merge_extr_stats_python} \
        --detected_tsv ${det_umi} \
        --synthetic_tsv ${extr_syn} \
        --umi_tsv ${extr_umi} \
        $write_report \
        $cons \
        -o .
    """
}