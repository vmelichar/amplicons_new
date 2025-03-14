process DETECT_UMI_CONSENSUS_FASTQ {
    // publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.output}/${sample}/${target}/${params.output_format}_umi/${type}", pattern: "*${params.output_format}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( fastq )
        val ( type )
        path umi_extract_python
    
    output:
        tuple val( "${sample}" ), val( "${target}" ), path ( "*${params.output_format}" ), emit: umi_extract_fastq
        tuple val( "${sample}" ), val( "${target}" ), path ( "*_detected_umis.tsv" ), path ( "*_extr_synthetic.tsv" ), path ( "*_extr_umi.tsv" ), emit: stats_tsv

    script:
        def write_report = params.write_reports ? "--tsv" : ""
        def cons = "${type}" == "consensus" ? "--cons" : ""

    """
        python ${umi_extract_python} \
        --fwd-umi ${params.fwd_umi} \
        --rev-umi ${params.rev_umi} \
        --max-error ${params.umi_errors} \
        --adapter_length ${params.adapter_length} \
        --output_format ${params.output_format} \
        --output_filename ${fastq.baseName}_detected_umis \
        --output_synthetic ${fastq.baseName}_extr_synthetic \
        --output_umi ${fastq.baseName}_extr_umi \
        $write_report \
        $cons \
        -o . ${fastq}
    """
}