process FILTER_CLUSTERS_PARALLEL {
    tag "${barcode}-${target}-batch${batch_id}"
    label 'process_low'
    cpus 2
    
    input:
        tuple val(barcode), val(target), val(batch_id), path(clusters)
        val(min_reads_per_cluster)
    
    output:
        tuple val(barcode), val(target), path("filtered/*.fasta"), optional: true, emit: filtered
        tuple val(barcode), val(target), val(low_count), emit: low_count
    
    script:
    """
    mkdir -p filtered
    LOW_COUNT=0
    
    # Parallel processing with GNU parallel if available
    if command -v parallel &> /dev/null; then
        export MIN_READS=${min_reads_per_cluster}
        
        parallel -j ${task.cpus} '
            count=\$(grep -c "^>" {} 2>/dev/null || echo 0)
            if [ \$count -ge \$MIN_READS ]; then
                cp {} filtered/
            else
                echo \$count
            fi
        ' ::: ${clusters} | awk '{sum+=\$1} END {print sum}' > low_count.txt
        
        LOW_COUNT=\$(cat low_count.txt)
    else
        # Fallback to sequential
        for cluster in ${clusters}; do
            count=\$(grep -c "^>" "\$cluster" 2>/dev/null || echo 0)
            if [ \$count -ge ${min_reads_per_cluster} ]; then
                cp "\$cluster" filtered/
            else
                LOW_COUNT=\$((LOW_COUNT + count))
            fi
        done
    fi
    
    echo \$LOW_COUNT > low_count.txt
    """
    
    stub:
    """
    mkdir -p filtered
    echo "0" > low_count.txt
    """
}