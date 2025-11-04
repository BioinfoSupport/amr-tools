
include { ORGANIZE_FILES } from './modules/organize_files'
include { NANOPLOT       } from './modules/nanoplot'
include { FASTQC         } from './modules/fastqc'
include { FASTP          } from './modules/fastp'
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

		// Run fastp once on each FASTQ-pair, and then expand to inputs sharing the same FASTQ
		FASTP(fqs_ch.map({m,x -> [x,x]}).unique())
		fastp_json_ch = fqs_ch.combine(FASTP.out.json,by:[1,0]).map({m,x,y -> [m,y]})
		fastp_html_ch = fqs_ch.combine(FASTP.out.html,by:[1,0]).map({m,x,y -> [m,y]})
		
		// Run FASTQC
		// WARNING: similarly to above, it should run once on each FASTQ and then expand to all inputs
		//          but it is not possible as we need to give sample_id to have correct sample names for MULTIQC
		//FASTQC(fqs_ch.map({meta,x -> [meta,x,meta.sample_id]}))

		// MultiQC
		ORGANIZE_FILES(
			Channel.empty().mix(
				nanostat_ch.map({meta,file -> [file,"${meta.sample_id}.nanostat"]}),
				//FASTQC.out.zip.map({meta,file -> [file,"${meta.sample_id}_fastqc.zip"]}),
				fastp_json_ch.map({meta,file -> [file,"${meta.sample_id}.fastp.json"]})
			)
			.collect({[it]})
			.map({["multiqc.html",it]})
		)
		MULTIQC(ORGANIZE_FILES.out,Channel.value(file("${moduleDir}/assets/multiqc_config.yml")))
	emit:
		long_nanoplot     = nanoplot_ch
		long_nanostat     = nanostat_ch
		//short_fastqc_html = FASTQC.out.html
		//short_fastqc_zip  = FASTQC.out.zip
		short_fastp_json  = fastp_json_ch
		short_fastp_html  = fastp_html_ch
		multiqc_html      = Channel.empty() //MULTIQC.out.html.map({m,x -> x})
}


workflow {
	main:
		println("TEST")
	publish:
		seq_multiqc = Channel.empty() //SEQUENCING_QC.out.multiqc_html
}

output {
	seq_multiqc {
		path { it >> "./multiqc_sequencing.html" }
	}
}
