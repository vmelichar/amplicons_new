process STRAND_STATS {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( strand_umi, list: true ), path ( strand_conca, list: true ), path ( strand_short, list: true ), path ( strand_long, list: true ), path ( strand_filter, list: true )
        val ( type )
    
    output:
        tuple val( "${sample}" ), val( "${target}" ), path ( "*.txt" )

    script:

    """
        grep strand=+ ${strand_conca[@]} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_conca[@]} | wc -l >> "strand_filter.txt"

        grep strand=+ ${strand_short[@]} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_short[@]} | wc -l >> "strand_filter.txt"

        grep strand=+ ${strand_long[@]} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_long[@]} | wc -l >> "strand_filter.txt"

        grep strand=+ ${strand_filter[@]} | wc -l >> "strand_filter.txt"
        grep strand=- ${strand_filter[@]} | wc -l >> "strand_filter.txt"
        
        grep strand=+ ${strand_umi[@]} | wc -l >> "umi_strand_filter.txt"
        grep strand=- ${strand_umi[@]} | wc -l >> "umi_strand_filter.txt"

    """
}