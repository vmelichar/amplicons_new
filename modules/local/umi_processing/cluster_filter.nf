process FILTER_CLUSTERS_PARALLEL {
    tag "${barcode}_${target}_batch${batch_idx}"
    label 'process_low'
    cpus 2
    input:
    tuple val(barcode), val(target), val(batch_idx), path(clusters)
    val(min_reads)
    
    output:
    tuple val(barcode), val(target), path("filtered/*.fasta"), optional: true, emit: filtered
    tuple val(barcode), val(target), env(LOW_COUNT), emit: low_count
    
    script:
    """
    #!/bin/bash -ue
    mkdir -p filtered
    export MIN_READS=${min_reads}
    
    # Create a file list for parallel
    ls *.fasta > cluster_files.txt 2>/dev/null || touch cluster_files.txt
    
    if command -v parallel &> /dev/null && [ -s cluster_files.txt ]; then
        # Use GNU parallel with proper quoting
        LOW_COUNT=\$(parallel -j ${task.cpus} '
            count=\$(grep -c "^>" {} 2>/dev/null || echo 0)
            if [ \$count -ge \$MIN_READS ]; then
                cp {} filtered/
                echo 0
            else
                echo \$count
            fi
        ' :::: cluster_files.txt | awk '{sum+=\$1} END {print sum+0}')
    else
        # Fallback to sequential
        LOW_COUNT=0
        for cluster in *.fasta; do
            [ -f "\$cluster" ] || continue
            count=\$(grep -c "^>" "\$cluster" 2>/dev/null || echo 0)
            if [ \$count -ge ${min_reads} ]; then
                cp "\$cluster" filtered/
            else
                LOW_COUNT=\$((LOW_COUNT + count))
            fi
        done
    fi
    
    echo "Batch ${batch_idx} filtered. Low-count reads: \$LOW_COUNT" >&2
    """
}