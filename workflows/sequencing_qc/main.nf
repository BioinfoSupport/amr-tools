
include { ORGANIZE_FILES } from './modules/organize_files'
include { NANOPLOT       } from './modules/nanoplot'
include { FASTQC         } from './modules/fastqc'
include { MULTIQC        } from './modules/multiqc'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow SEQUENCING_QC {
	take:
		fql_ch    // channel: [ val(meta), path(long_reads) ]
		fqs_ch    // channel: [ val(meta), path(short_reads) ]
	main:
		// Reads Quality Controls
		
		// Run nanoplot once on each FASTQ, and then expand to inputs sharing the same FASTQ
		NANOPLOT(fql_ch.map({m,x -> [x,x]}).unique())
		nanostat_ch = fql_ch.combine(NANOPLOT.out.nanostat,by:[1,0]).map({m,x,y -> [m,y]})
		nanoplot_ch = fql_ch.combine(NANOPLOT.out.nanoplot,by:[1,0]).map({m,x,y -> [m,y]})
		
		// Run FASTQC
		// WARNING: similarly to above, it should run once on each FASTQ and then expand to all inputs
		//          but it is not possible as we need to give sample_id to have correct sample names for MULTIQC
		FASTQC(fqs_ch.map({meta,x -> [meta,x,meta.sample_id]}))

		// MultiQC
		ORGANIZE_FILES(
			Channel.empty().mix(
				nanostat_ch.map({meta,file -> [file,"${meta.sample_id}.nanostat"]}),
				FASTQC.out.zip.map({meta,file -> [file,"${meta.sample_id}_fastqc.zip"]})
			)
			.collect({[it]})
			.map({["multiqc.html",it]})
		)
		MULTIQC(ORGANIZE_FILES.out,file("${moduleDir}/assets/multiqc_config.yml"))

	emit:
		long_nanoplot     = nanoplot_ch
		long_nanostat     = nanostat_ch
		short_fastqc_html = FASTQC.out.html
		short_fastqc_zip  = FASTQC.out.zip
		multiqc_html      = MULTIQC.out.html.map({m,x -> x})
}


