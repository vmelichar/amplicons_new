include {MERGE_FASTQ} from '../modules/local/umi_processing/merge_input.nf'
include {SUBSAMPLING} from '../modules/local/umi_processing/subsampling.nf'
include {MAP_READS} from '../modules/local/umi_processing/map_reads.nf'
include {SPLIT_READS} from  '../modules/local/umi_processing/split_reads.nf'
include {DETECT_UMI_CONSENSUS_FASTQ as DETECT_UMI_FASTQ} from '../modules/local/umi_polishing/detect_umi_consensus_fastq.nf'
include {CLUSTER} from '../modules/local/umi_processing/cluster.nf'
include {REFORMAT_FILTER_CLUSTER} from '../modules/local/umi_processing/reformat_filter_cluster.nf'
include {CLUSTER_STATS} from '../modules/local/umi_processing/cluster_stats.nf'
include {SUMMARY_CLUSTER_STATS} from '../modules/local/umi_processing/summary_cluster_stats.nf'
include {MERGE_EXTRACTION_STATS} from '../modules/local/umi_processing/merge_extraction_stats.nf'
include {MERGE_FILTER_STATS} from '../modules/local/umi_processing/merge_filter_stats.nf'
include {STRAND_STATS} from '../modules/local/umi_processing/strand_stats.nf'


workflow OFFLINE_UMI_PROCESSING {

    take:
        raw
        reference
        umi_filter_reads
        umi_extract
        umi_parse_clusters
        umi_cluster_report
        umi_cluster_stats_summary
        cluster_summary_cache_dir_nf

        bed
        merge_extr_stats
        merge_filter_stats

    main:       
        Channel
            .fromPath("${params.input}/barcode*/*.fastq")
            .map{ 
                fastqs -> 
                def barcode = fastqs.parent.name
                tuple(barcode, fastqs)
                }
            .groupTuple( by: 0 ) 
            .set{ existing_fastqs }

        if( params.subsampling ){
            MERGE_FASTQ( existing_fastqs )
            SUBSAMPLING( MERGE_FASTQ.out.merged_fastq, raw )
            SUBSAMPLING.out.subsampled_fastq
            .set { input_fastqs }
        } else {
            MERGE_FASTQ( existing_fastqs )
            .set { input_fastqs }
        }

        input_fastqs
            .splitFastq( by: params.chunk_size , file: true)
            .set{ chunked_input_fastqs }

        MAP_READS( chunked_input_fastqs, raw, reference )
        
        MAP_READS.out.bam_consensus
        .combine(bed)
        .set{ bam_consensus_bed_sets }

        SPLIT_READS( bam_consensus_bed_sets, raw, umi_filter_reads )

        SPLIT_READS.out.stats_tsv
        .groupTuple( by: [0, 1])
        .set{ filter_stats_to_merge }

        MERGE_FILTER_STATS( filter_stats_to_merge, raw, merge_filter_stats )

        SPLIT_READS.out.split_reads_fastx
        .filter{ _sample, _target, fastq -> fastq.countFastq() > params.min_reads_per_barcode }
        .set{ split_reads_filtered }

        DETECT_UMI_FASTQ( split_reads_filtered, raw, umi_extract )
        
        DETECT_UMI_FASTQ.out.umi_extract_fastq
        .groupTuple( by: [0, 1])
        .set{ extracted_umis }

        SPLIT_READS.out.split_reads_fastx_conca
        .groupTuple( by: [0, 1])
        .set{ strand_conca }
        
        SPLIT_READS.out.split_reads_fastx_short
        .groupTuple( by: [0, 1])
        .set{ strand_short }
        
        SPLIT_READS.out.split_reads_fastx_long
        .groupTuple( by: [0, 1])
        .set{ strand_long }

        SPLIT_READS.out.split_reads_fastx
        .groupTuple( by: [0, 1])
        .set{ strand_filter }

        extracted_umis
        .join(strand_conca, by: [0, 1])
        .join(strand_short, by: [0, 1])
        .join(strand_long, by: [0, 1])
        .join(strand_filter, by: [0, 1])
        .set { channel_to_strand }

        STRAND_STATS(channel_to_strand, raw)

        DETECT_UMI_FASTQ.out.stats_tsv
        .groupTuple( by: [0, 1])
        .set{ stats_to_merge }

        MERGE_EXTRACTION_STATS( stats_to_merge, raw, merge_extr_stats )

        CLUSTER( extracted_umis, raw )

        // Filters the clusters to only keep cluser with more or equal than min_reads_per_cluster, but keeps the grouping per sample
        CLUSTER.out.cluster_fastas
            .map { barcode, target, clusters -> 
                def filtered_clusters = clusters.findAll { fasta -> fasta.countFasta() >= params.min_reads_per_cluster }
                filtered_clusters ? [barcode, target, filtered_clusters] : null
            }
            .filter { it != null }
            .set{ cluster_fastas }

        CLUSTER.out.cluster_fastas
            .map { barcode, target, clusters -> 
                def total_low_count = clusters.findAll { fasta -> fasta.countFasta() < params.min_reads_per_cluster }
                                      .sum { fasta -> fasta.countFasta() } ?: 0
                total_low_count > 0 ? [barcode, target, total_low_count] : null
            }
            .filter { it != null }
            .groupBy { it[0] } // Group by barcode
            .map { barcode, values -> 
                def target_counts = values.collectEntries { [it[1], it[2]] } // {target: count}
                [barcode, target_counts] // Single unpackable value
                }
            .set { low_clusters_counts }

        REFORMAT_FILTER_CLUSTER( cluster_fastas, raw, umi_parse_clusters )

        CLUSTER_STATS( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, raw, umi_cluster_report )
        SUMMARY_CLUSTER_STATS( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, cluster_summary_cache_dir_nf, umi_cluster_stats_summary)

        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
        .filter{ _sample, _type, fastqs, _task_index -> fastqs instanceof List}            
        .map{ sample, type, fastqs, _task_index ->
                tuple(sample, type, fastqs)
        }
        .set{ processed_umis }

        emit:
            processed_umis
            low_clusters_counts

}