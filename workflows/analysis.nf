include {ANALYSIS_SNP} from '../modules/local/analysis/analysis_snp.nf'

workflow ANALYSIS {
    take:
        final_bam
        analysis_script
        positions_file

    main:
        ANALYSIS_SNP(final_bam, positions_file, analysis_script)
        

}