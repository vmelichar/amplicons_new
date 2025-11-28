process STRAND_STATS {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val(sample), val(target), path(strand_umi), path(strand_short), path(strand_long), path(strand_filter)
        val(type)

    output:
        tuple val(sample), val(target), path("*.txt")

    script:
    """
        # helper function to safely count strands
        count_strands() {
            files=\$1
            if [ -z "\$files" ]; then
                echo 0
                echo 0
            else
                grep -h 'strand=+' \$files | wc -l
                grep -h 'strand=-' \$files | wc -l
            fi
        }

        {
            count_strands "${strand_short.join(' ')}"
            count_strands "${strand_long.join(' ')}"
            count_strands "${strand_filter.join(' ')}"
        } > strand_filter.txt

        {
            count_strands "${strand_umi.join(' ')}"
        } > umi_strand_filter.txt
    """
}
