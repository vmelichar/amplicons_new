process RENAME_SEQUENCES {
    tag "${barcode}"
    
    input:
    tuple val(barcode), path(fastqs)
    
    output:
    tuple val(barcode), path("r_*.fastq.gz")
    
    script:
    def input_files = fastqs instanceof List ? fastqs : [fastqs]
    def commands = input_files.collect { fq ->
        def base_name = fq.name.replaceAll(/\.fastq(\.gz)?$/, '')
        def output_name = "r_${fq.name}"

        if (fq.name.endsWith('.gz')) {
            """
            zcat ${fq} | awk -v prefix="${base_name}" '
                NR % 4 == 1 {
                    if (\$0 ~ /^@/) {
                        hash = substr(\$0, 2)
                        print "@" prefix ":" hash
                    } else {
                        print \$0
                    }
                    next
                }
                { print }
            ' | gzip > ${output_name}
            """
        } else {
            """
            awk -v prefix="${base_name}" '
                NR % 4 == 1 {
                    if (\$0 ~ /^@/) {
                        hash = substr(\$0, 2)
                        print "@" prefix ":" hash
                    } else {
                        print \$0
                    }
                    next
                }
                { print }
            ' ${fq} > ${output_name}
            """
        }
    }.join('\n    ')
    
    """
    ${commands}
    """
}