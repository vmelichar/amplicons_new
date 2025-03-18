process STRAND_STATS {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( strand_umi )
        tuple val( sample ), val( target ), path ( strand_conca )
        tuple val( sample ), val( target ), path ( strand_short )
        tuple val( sample ), val( target ), path ( strand_long )
        tuple val( sample ), val( target ), path ( strand_filter )
        val ( type )
    
    output:
        tuple val( "${sample}" ), val( "${target}" ), path ( "*.txt" )

    script:

    """
        grep strand=+ ${strand_conca} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_conca} | wc -l >> "strand_filter.txt"

        grep strand=+ ${strand_short} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_short} | wc -l >> "strand_filter.txt"

        grep strand=+ ${strand_long} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_long} | wc -l >> "strand_filter.txt"

        grep strand=+ ${strand_filter} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_filter} | wc -l >> "strand_filter.txt"
        
        grep strand=+ ${strand_umi} | wc -l >> "umi_strand_filter.txt"
        grep strand=- ${strand_umi} | wc -l >> "umi_strand_filter.txt"

    """
}