process MERGE_EXTRACTION_STATS {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "merged_det_umi.tsv", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( det_umi ), path ( extr_syn ), path ( extr_umi )
        val ( type )
        path merge_extr_stats_python
    
    output:
        path "merged_det_umi.tsv"

    script:
        def write_report = params.write_reports ? "--tsv" : ""
        def cons = "${type}" == "consensus" ? "--cons" : ""

    """
        python ${merge_extr_stats_python} \
        --detected_tsv ${det_umi} \
        --synthetic_tsv ${extr_syn} \
        --umi_tsv ${extr_umi} \
        $write_report
        $cons
        -o .
    """
}