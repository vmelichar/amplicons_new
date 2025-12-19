process MOLECULE_STATS {
    tag "${sample}_${target}"
    publishDir "${params.output}/${sample}/${target}/stats/", pattern: "*.tsv", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( extr_syn ), path ( extr_umi ), path ( cluster_hash ), path ( extr_syn_cons ), path ( extr_umi_cons ), path ( seq_type_hash )
        path molecule_stats_python
    
    output:
        path "*_molecule_statistics.tsv", emit: stats_tsv

    script:
    """
        python ${molecule_stats_python} \
        --hotspot ${target} \
        --extraction_raw_syn_file ${extr_syn} \
        --extraction_raw_umi_file ${extr_umi} \
        --cluster_hash_file ${cluster_hash} \
        --seq_type_hash_file ${seq_type_hash} \
        --extraction_cons_syn_file ${extr_syn_cons} \
        --extraction_cons_umi_file ${extr_umi_cons} \
        -o .
    """
}