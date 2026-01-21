
include { SAMTOOLS_FASTQ } from './modules/samtools/fastq'
include { NANOPLOT       } from './modules/nanoplot'
include { FASTP          } from './modules/fastp'
include { ORGANIZE_FILES } from './modules/organize_files'
include { MULTIQC        } from './modules/multiqc'



workflow SEQ_QC {
	take:
		fqs_ch    // channel: [ val(meta), path(short_reads) ]	
		fql_ch    // channel: [ val(meta), path(long_reads) ]
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
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

/*
process GENERATE_READSETS {
  output:
    path "readsets.csv"
  script:
  	//def lr = Channel.fromPath(long_reads_pattern)
	  """
	  touch readsets.csv
	  """
}
*/


params.readsets = null
params.long_reads = null
params.short_reads = null

import AmrUtils

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))
/*
		if (params.readsets) {
				readsets = Channel.fromPath(params.readsets,checkIfExists:true)
		} else {
		    readsets = GENERATE_READSETS()
		}
		println(readsets.collect())
		//println(readsets.map({samplesheetToList(it, "assets/schema_readsets.json")}))
		//Channel.fromList(samplesheetToList(readsets, "assets/schema_readsets.json")).view()
*/		
		
		def ss = AmrUtils.parse_generic_params(params,{sheet -> samplesheetToList(sheet, "assets/schema_samplesheet.json")})

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		ss.fql_ch = ss.fql_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		ss.fql_ch = ss.fql_ch.fq.mix(SAMTOOLS_FASTQ(ss.fql_ch.bam))

		// Reduce FASTQ size if needed
		//lr_ch = FQ_SUBSAMPLE(ss.lr_ch)

		// Reads Quality Controls, get a multiQC
		SEQ_QC(ss.fqs_ch,ss.fql_ch)
		
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





