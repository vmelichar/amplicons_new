process STRAND_STATS {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( strand_umi ), path ( strand_conca ), path ( strand_short ), path ( strand_long ), path ( strand_filter )
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