include {GLUE_CLUSTERS} from '../modules/local/umi_polishing/glue_clusters.nf'
include {POLISH_CLUSTER} from '../modules/local/umi_polishing/polish_cluster.nf'
include {MERGE_CONSENSUS_FASTQ} from '../modules/local/umi_polishing/merge_consensus_fastq.nf'
include {FILTER_CONSENSUS_FASTQ} from '../modules/local/umi_polishing/filter_consensus_fastq.nf'
include {REFORMAT_CONSENSUS_CLUSTER} from '../modules/local/umi_polishing/reformat_consensus_cluster.nf'
include {MAP_CONSENSUS; MAP_CONSENSUS as MAP_FINAL_CONSENSUS} from '../modules/local/umi_polishing/map_consensus.nf'
include {DETECT_UMI_CONSENSUS_FASTQ} from '../modules/local/umi_polishing/detect_umi_consensus_fastq.nf'
include {CLUSTER_CONSENSUS} from '../modules/local/umi_polishing/cluster_consensus.nf'
include {MERGE_EXTRACTION_STATS as MERGE_CONSENSUS_EXTRACTION_STATS} from '../modules/local/umi_processing/merge_extraction_stats.nf'

workflow UMI_POLISHING {
    take:
        processed_umis
        n_parsed_cluster
        consensus
        final_consensus
        reference
        umi_extract
        umi_reformat_consensus
        merge_extr_stats
        bed_ch

    main:
         GLUE_CLUSTERS(processed_umis, consensus)
            .map{ sample, target, clusters -> tuple(sample, target, clusters instanceof List ? clusters : [clusters]) }
            .set{ glued_clusters }
        
        glued_clusters
        .map{ sample, _target, clusters -> n_parsed_cluster.put("$sample", clusters.size)}
        
        glued_clusters
            .transpose(by: 2)
            .set { glued_clusters_transposed }

        POLISH_CLUSTER( glued_clusters_transposed, consensus )
        
        POLISH_CLUSTER.out.consensus_fastq
        .map{ sample, target, fastq -> tuple( groupKey(sample, n_parsed_cluster.get("$sample")), target, fastq) }
        .groupTuple( by: [0,1] )
        .set{ merge_consensus }

        
        if ( params.output_format == "fastq"){
            MERGE_CONSENSUS_FASTQ(merge_consensus, consensus)
            FILTER_CONSENSUS_FASTQ(MERGE_CONSENSUS_FASTQ.out.merged_consensus_fastq, consensus)
            FILTER_CONSENSUS_FASTQ.out.filtered_consensus_fastq
            .set{ consensus_fastq }
        } else {
            MERGE_CONSENSUS_FASTQ(merge_consensus, consensus)
            .set{ consensus_fastq }
        }

        MAP_CONSENSUS( consensus_fastq, consensus, reference )
        DETECT_UMI_CONSENSUS_FASTQ( consensus_fastq, consensus, umi_extract )

        DETECT_UMI_CONSENSUS_FASTQ.out.umi_extract_fastq
            .map { barcode, target, fastq -> 
                def absoluteFastq = fastq.collect { it.toAbsolutePath() }
                def filtered_fastq = absoluteFastq.findAll { fastq -> fastq.countFastq() > 0 }
                filtered_fastq ? [barcode, target, filtered_fastq] : null
            }
            .filter { it != null }
            .set{ detected_fastq }

        DETECT_UMI_CONSENSUS_FASTQ.out.stats_tsv
        .groupTuple( by: [0, 1])
        .set{ stats_to_merge_cons }

        MERGE_CONSENSUS_EXTRACTION_STATS( stats_to_merge_cons, consensus, merge_extr_stats )

        CLUSTER_CONSENSUS( detected_fastq , consensus )
        REFORMAT_CONSENSUS_CLUSTER( CLUSTER_CONSENSUS.out.consensus_fasta, final_consensus, umi_reformat_consensus )
        MAP_FINAL_CONSENSUS( REFORMAT_CONSENSUS_CLUSTER.out.consensus_fastq, final_consensus, reference )
        
        MAP_FINAL_CONSENSUS.out.bam_consensus
        .map{ sample, target, bam, bai -> tuple(target, sample, bam, bai)}
        .join(bed_ch)
        .set{snp_analysis_bam}

    emit: 
        consensus_bam = MAP_CONSENSUS.out.bam_consensus
        final_consensus_bam = MAP_FINAL_CONSENSUS.out.bam_consensus
        snp_analysis_bam

}
