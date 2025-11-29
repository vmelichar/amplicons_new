include {ANALYSIS_SNP} from '../modules/local/analysis/analysis_snp.nf'
include {EXPORT_FILTER_STATS} from '../modules/local/umi_polishing/export_filter_stats.nf'

workflow ANALYSIS {
    take:
        final_bam
        analysis_script
        positions_file
        variants_vcf
        variants_tbi
        export_stats
        low_clusters_counts

    main:
        ANALYSIS_SNP(final_bam, positions_file, variants_vcf, variants_tbi, analysis_script)
        
        EXPORT_FILTER_STATS( export_stats, low_clusters_counts, ANALYSIS_SNP.out.flag_file )

}