
include { SAMTOOLS_FASTQ    } from './modules/samtools/fastq'
include { NANOPLOT       } from './modules/nanoplot'
include { FASTP          } from './modules/fastp'
include { ORGANIZE_FILES } from './modules/organize_files'
include { MULTIQC        } from './modules/multiqc'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'



workflow SEQ_QC {
	take:
		fql_ch    // channel: [ val(meta), path(long_reads) ]
		fqs_ch    // channel: [ val(meta), path(short_reads) ]
	main:
		// Run nanoplot once on each FASTQ, and then expand to inputs sharing the same FASTQ
		NANOPLOT(fql_ch)

		// Run fastp once on each FASTQ-pair, and then expand to inputs sharing the same FASTQ
		FASTP(fqs_ch)

		// MultiQC
		ORGANIZE_FILES(
			Channel.empty().mix(
				NANOPLOT.out.nanostat.map({meta,file -> [file,"${meta.sample_id}.nanostat"]}),
				//FASTQC.out.zip.map({meta,file -> [file,"${meta.sample_id}_fastqc.zip"]}),
				FASTP.out.json.map({meta,file -> [file,"${meta.sample_id}.fastp.json"]})
			)
			.collect({[it]})
			.map({["multiqc.html",it]})
		)
		MULTIQC(ORGANIZE_FILES.out,Channel.value(file("${moduleDir}/assets/multiqc_config.yml")))
	emit:
		long_nanoplot     = NANOPLOT.out.nanoplot
		long_nanostat     = NANOPLOT.out.nanostat
		short_fastp_json  = FASTP.out.json
		short_fastp_html  = FASTP.out.html
		multiqc_html      = MULTIQC.out.html.map({m,x -> x})
}





// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------

params.samplesheet = null
params.long_reads = null
params.short_reads = null

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Extract Long_read and Short_read channels from params 
		lr_ch = Channel.empty()
		sr_ch = Channel.empty()
		if (params.samplesheet) {
			SS = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
				.multiMap({x ->
					lr_ch: [[sample_id:x[0].sample_id],x[0].long_reads]
					sr_ch: [[sample_id:x[0].sample_id],[x[0].short_reads_1,x[0].short_reads_2]]
				})
			lr_ch = SS.lr_ch
			sr_ch = SS.sr_ch
		} else {
			if (params.long_reads) {
				lr_ch = Channel.fromPath(params.long_reads)
						.map({x -> tuple(["sample_id":x.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')],x)})
			}
			if (params.short_reads) {
				sr_ch = Channel
						.fromFilePairs(params.short_reads,size:-1) { file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') }
						.map({id,x -> [["sample_id":id],x]})
			}
		}
		// Filter out missing values
		sr_ch = sr_ch.map({x,y -> [x,y.findAll({v->v})]}).filter({x,y -> y})
		lr_ch = lr_ch.filter({x,y -> y})

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		lr_ch = lr_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		lr_ch = lr_ch.fq.mix(SAMTOOLS_FASTQ(lr_ch.bam))

		// Reduce FASTQ size if needed
		//lr_ch = FQ_SUBSAMPLE(ss.lr_ch)

		// Reads Quality Controls, get a multiQC
		SEQ_QC(lr_ch,sr_ch)
		
	publish:
		long_nanoplot     = SEQ_QC.out.long_nanoplot
		short_fastp_json  = SEQ_QC.out.short_fastp_json
		short_fastp_html  = SEQ_QC.out.short_fastp_html
		multiqc_html      = SEQ_QC.out.multiqc_html
}

output {
	long_nanoplot {
		path { m,x -> x >> "samples/${m.sample_id}/seq_qc/long_nanoplot"}
	}
	short_fastp_json {
		path { m,x -> x >> "samples/${m.sample_id}/seq_qc/short_fastp.json"}
	}	
	short_fastp_html {
		path { m,x -> x >> "samples/${m.sample_id}/seq_qc/short_fastp.html"}
	}
	multiqc_html {
		path { it >> "./seq_qc.html" }
	}
}





