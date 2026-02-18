

include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_LONG  } from './modules/fq_subsample'
include { CHOPPER as CHOPPER_LONG            } from './modules/chopper'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_SHORT } from './modules/fq_subsample'


workflow READS_FILTER {
	take:
		fqs_ch    // channel: [ val(meta), path(short_reads) ]	
		fql_ch    // channel: [ val(meta), path(long_reads) ]
	main:
		// Reduce FASTQ size if needed
		fql_ch = fql_ch | CHOPPER_LONG | FQ_SUBSAMPLE_LONG
		fqs_ch = fqs_ch | FQ_SUBSAMPLE_SHORT
	emit:
	  long_filtered     = fql_ch
		short_filtered    = fqs_ch
}

