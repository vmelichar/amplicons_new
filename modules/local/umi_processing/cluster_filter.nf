process FILTER_CLUSTERS_PARALLEL {
    tag "${barcode}_${target}_batch${batch_idx}"
    cpus 2
    
    input:
    tuple val(barcode), val(target), val(batch_idx), path('cluster_input/*')  // Stage all paths in the list to this directory
    
    output:
    tuple val(barcode), val(target), path("filtered/*.fasta"), optional: true, emit: filtered
    tuple val(barcode), val(target), env(LOW_COUNT), emit: low_count
    
    script:
    """
    #!/bin/bash -ue
    mkdir -p filtered
    export MIN_READS=${params.min_reads_per_cluster}

    echo "Processing batch ${batch_idx} for ${barcode}/${target}"
    FILE_COUNT=\$(ls cluster_input/ | wc -l)
    echo "Files in batch: \$FILE_COUNT"

    if command -v parallel &> /dev/null; then
        # Use GNU parallel
        LOW_COUNT=\$(ls cluster_input/ | parallel -j ${task.cpus} '
            cluster="cluster_input/{}"
            count=\$(grep -c "^>" "\$cluster" 2>/dev/null || echo 0)
            if [ \$count -ge '"${params.min_reads_per_cluster}"' ]; then
                cp "\$cluster" filtered/
                echo 0
            else
                echo \$count
            fi
        ' | awk '{sum+=\$1} END {print sum+0}')
    else
        # Fallback sequential
        LOW_COUNT=0
        for cluster in cluster_input/*; do
            [ -f "\$cluster" ] || continue
            count=\$(grep -c "^>" "\$cluster" 2>/dev/null || echo 0)
            if [ \$count -ge ${params.min_reads_per_cluster} ]; then
                cp "\$cluster" filtered/
            else
                LOW_COUNT=\$((LOW_COUNT + count))
            fi
        done
    fi

    echo "Batch ${batch_idx} complete: \$(ls filtered/ 2>/dev/null | wc -l) passed, \$LOW_COUNT low-count reads"
    """
}