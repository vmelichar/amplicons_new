process ANALYSIS_SNP {
    tag "${sample}_${target}"
    conda '/home/melichv/miniconda3/envs/amplicons'
    publishDir "${params.output}/${sample}/${target}/analysis", pattern: "*.csv", mode: 'copy'
    publishDir "${params.output}/${sample}/${target}/analysis", pattern: "*.png", mode: 'copy'

    input:
        tuple val( target ), val( sample ), path ( bam ), path ( bai ), path ( bed )
        path positions_file
        path variants_vcf
        path variants_tbi
        path analysis_script
    
    output:
        path "*.csv"
        path "*.png"
        path "*.txt"
        tuple val( "${sample}" ), val( "${target}" ), path( "*_cluster_types.tsv" ), emit: cluster_types_file
        path "flag_file.md", emit: flag_file

    script:

    """
        python ${analysis_script} \
        --hs ${target} \
        --bam ${bam} \
        --bai ${bai} \
        --vcf ${variants_vcf} \
        --tbi ${variants_tbi} \
        --positions ${positions_file} \
        -o .
        
        touch flag_file.md
    """
}