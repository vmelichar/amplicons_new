include {ANALYSIS_SNP} from '../modules/local/analysis/analysis_snp.nf'

workflow ANALYSIS {
    take:
        final_bam
        analysis_script
        positions_file
        variants_vcf
        variants_tbi

    main:
        ANALYSIS_SNP(final_bam, positions_file, variants_vcf, variants_tbi, analysis_script)
        

}