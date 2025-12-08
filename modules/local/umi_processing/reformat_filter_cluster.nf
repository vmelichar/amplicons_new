process REFORMAT_FILTER_CLUSTER {
    beforeScript "rm -f ${params.output}/${sample}/${target}/clustering/${type}/*.tar.gz"
    tag "${sample}_${target}"
    publishDir "${params.output}/${sample}/${target}/clustering/${type}", pattern: "*.tar.gz", mode: 'copy'
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*tsv", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path( cluster )
        val( type )
        path umi_parse_clusters_python

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "smolecule*"), val(task.index), optional: true, emit: smolecule_cluster_fastqs
        tuple val( "${sample}" ), val ( "${target}" ), path( "*.tsv" ), optional: true, emit: smolecule_cluster_stats
        tuple val( "${sample}" ), val( "${target}" ), path( "*.tar.gz"), optional: true

    script:
        def balance_strands = params.balance_strands ? "--balance_strands" : ""
        def write_report = params.write_reports ? "--tsv" : ""
    """
        # Save input cluster file paths to a temporary list file
        echo ${cluster} | tr ' ' '\n' > cluster_list.txt

        python ${umi_parse_clusters_python} \
          --filter_strategy ${params.filter_strategy_clusters} \
          --min_reads_per_clusters ${params.min_reads_per_cluster} \
          --max_reads_per_clusters ${params.max_reads_per_cluster} \
          --max_dist_umi ${params.max_dist_umi} \
          --cluster_list cluster_list.txt \
          --output_format ${params.output_format} \
          $balance_strands \
          $write_report \
          -o .

        tar -czf UMI_smolecules.tar.gz UMI_smolecule* || true
    """
}
