
include { SAMTOOLS_FASTQ     } from './modules/samtools/fastq'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_LONG  } from './modules/fq_subsample'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_SHORT } from './modules/fq_subsample'
include { NANOPLOT           } from './modules/nanoplot'
include { FASTP              } from './modules/fastp'
include { ORGANIZE_FILES     } from './modules/organize_files'
include { MULTIQC            } from './modules/multiqc'



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
		fql_ch = FQ_SUBSAMPLE_LONG(fql_ch,params.limit_long_reads_len)
		fqs_ch = FQ_SUBSAMPLE_SHORT(fqs_ch,params.limit_short_reads_len)

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
		short_filtered    = fqs_ch
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
params.limit_long_reads_len  = 1000000000
params.limit_short_reads_len = 1000000000

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))
		
		def readsets = Readsets.fromParams(params.readsets,{sheet,schema -> samplesheetToList(sheet, schema)})
		
		SEQ_QC(params,Readsets.short_reads_channel(readsets),Readsets.long_reads_channel(readsets))
  	
	publish:
	  input_readsets_csv    = Readsets.toCSV(readsets)
	  filtered_readsets_csv = Readsets.toCSV(Readsets.fromChannels(SEQ_QC.out.short_filtered,SEQ_QC.out.long_filtered))
		long_nanoplot         = SEQ_QC.out.long_nanoplot
		short_fastp_json      = SEQ_QC.out.short_fastp_json
		short_fastp_html      = SEQ_QC.out.short_fastp_html
		multiqc_html          = SEQ_QC.out.multiqc_html
}

output {
	input_readsets_csv {
    index {
    	path 'indexes/seq_qc_input_readsets.csv'
    	header true
    }
	}
	filtered_readsets_csv {
		path { m -> 
			m.long_reads >> "samples/${m.sample_id}/seq_qc/${m.readset_id}/filtered_long_reads.fastq.gz"
			m.short_reads_1 >> "samples/${m.sample_id}/seq_qc/${m.readset_id}/filtered_short_reads_1.fastq.gz"
			m.short_reads_2 >> "samples/${m.sample_id}/seq_qc/${m.readset_id}/filtered_short_reads_2.fastq.gz"
		}
    index {
    	path 'indexes/seq_qc_filtered_readsets.csv'
    	header true
    }
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





