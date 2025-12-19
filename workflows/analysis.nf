include {ANALYSIS_SNP} from '../modules/local/analysis/analysis_snp.nf'
include {EXPORT_FILTER_STATS} from '../modules/local/umi_polishing/export_filter_stats.nf'
include {MOLECULE_STATS} from '../modules/local/analysis/molecule_stats.nf'

workflow ANALYSIS {
    take:
        final_bam
        analysis_script
        positions_file
        variants_vcf
        variants_tbi
        export_stats
        low_clusters_counts
        extr_synthetic_stats
        extr_umi_stats
        cluster_read_hash
        extr_synthetic_stats_cons
        extr_umi_stats_cons
        molecule_stats_python

    main:
        ANALYSIS_SNP(final_bam, positions_file, variants_vcf, variants_tbi, analysis_script)
        
        EXPORT_FILTER_STATS( export_stats, low_clusters_counts, ANALYSIS_SNP.out.flag_file )

        extr_synthetic_stats
            .join( extr_umi_stats, by: [0,1] )
            .join( cluster_read_hash, by: [0,1] )
            .join( extr_synthetic_stats_cons, by: [0,1] )
            .join( extr_umi_stats_cons, by: [0,1] )
            .join( ANALYSIS_SNP.out.cluster_types_file, by: [0,1] )
            .set { molecule_stats_input }

        MOLECULE_STATS( molecule_stats_input, molecule_stats_python )
}