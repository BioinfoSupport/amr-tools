
include { SAMTOOLS_FASTQ     } from './modules/samtools/fastq'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_LONG  } from './modules/fq_subsample'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_SHORT } from './modules/fq_subsample'
include { NANOPLOT           } from './modules/nanoplot'
include { FASTP              } from './modules/fastp'
include { ORGANIZE_FILES     } from './modules/organize_files'
include { MULTIQC            } from './modules/multiqc'



workflow READS_FILTER {
	take:
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
		fql_ch = FQ_SUBSAMPLE_LONG(fql_ch)
		fqs_ch = FQ_SUBSAMPLE_SHORT(fqs_ch)

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




