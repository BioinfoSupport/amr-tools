
include { SAMTOOLS_FASTQ } from './modules/samtools/fastq'
include { FQ_SUBSAMPLE   } from './modules/fq_subsample'
include { NANOPLOT       } from './modules/nanoplot'
include { FASTP          } from './modules/fastp'
include { ORGANIZE_FILES } from './modules/organize_files'
include { MULTIQC        } from './modules/multiqc'



workflow SEQ_QC {
	take:
		opts
		fqs_ch    // channel: [ val(meta), path(short_reads) ]	
		fql_ch    // channel: [ val(meta), path(long_reads) ]
	main:
		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		fql_ch = fql_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		fql_ch = fql_ch.fq.mix(SAMTOOLS_FASTQ(fql_ch.bam))

		// Reduce FASTQ size if needed
		if (opts.limit_long_reads_len) {
			fql_ch = FQ_SUBSAMPLE(fql_ch,params.limit_long_reads_len)
		}

		// Run nanoplot once on each FASTQ, and then expand to inputs sharing the same FASTQ
		NANOPLOT(fql_ch)

		// Run fastp once on each FASTQ-pair, and then expand to inputs sharing the same FASTQ
		FASTP(fqs_ch)

		// MultiQC
		ORGANIZE_FILES(
			Channel.empty().mix(
				NANOPLOT.out.nanostat.map({meta,file -> [file,"${meta.sample_id}_${meta.readset_id}.nanostat"]}),
				//FASTQC.out.zip.map({meta,file -> [file,"${meta.sample_id}_${meta.readset_id}_fastqc.zip"]}),
				FASTP.out.json.map({meta,file -> [file,"${meta.sample_id}_${meta.readset_id}.fastp.json"]})
			)
			.collect({[it]})
			.map({["multiqc.html",it]})
		)
		MULTIQC(ORGANIZE_FILES.out,Channel.value(file("${moduleDir}/assets/multiqc_config.yml")))
	emit:
	  long_filtered     = fql_ch
		long_nanoplot     = NANOPLOT.out.nanoplot
		long_nanostat     = NANOPLOT.out.nanostat
		short_fastp_json  = FASTP.out.json
		short_fastp_html  = FASTP.out.html
		multiqc_html      = MULTIQC.out.html.map({m,x -> x})
}




// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
import Readsets
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.readsets = [
	csv         : null,
	long_reads  : [],
	short_reads : []
]
params.limit_long_reads_len = null

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		def readsets = new Readsets(params.readsets,{sheet,schema -> samplesheetToList(sheet, schema)})
		def input_readsets_csv = readsets.flat_csv()
		SEQ_QC(params,readsets.short_reads,readsets.long_reads)
    
	publish:
	  input_readsets_csv  = input_readsets_csv
	  filtered_long_reads = SEQ_QC.out.long_filtered
		long_nanoplot       = SEQ_QC.out.long_nanoplot
		short_fastp_json    = SEQ_QC.out.short_fastp_json
		short_fastp_html    = SEQ_QC.out.short_fastp_html
		multiqc_html        = SEQ_QC.out.multiqc_html
}

output {
	input_readsets_csv {
    index {
    	path 'indexes/seq_qc_input_readsets.csv'
    	header true
    }
	}
	filtered_long_reads {
		path { m,x -> x >> "samples/${m.sample_id}/seq_qc/${m.readset_id}/filtered_long_reads.fastq.gz"}
	}
	long_nanoplot {
		path { m,x -> x >> "samples/${m.sample_id}/seq_qc/${m.readset_id}/long_nanoplot"}
	}
	short_fastp_json {
		path { m,x -> x >> "samples/${m.sample_id}/seq_qc/${m.readset_id}/short_fastp.json"}
	}	
	short_fastp_html {
		path { m,x -> x >> "samples/${m.sample_id}/seq_qc/${m.readset_id}/short_fastp.html"}
	}
	multiqc_html {
		path { it >> "./seq_qc.html" }
	}
}





