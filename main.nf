

include {READS_QC}      from './wf/reads_qc/main.nf'
include {ASSEMBLE}      from './wf/assemble/main.nf'
include {ASSEMBLIES_QC} from './wf/assemblies_qc/main.nf'
include {AMR_ANNOT}     from './wf/amr_annot/main.nf'


// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
import Samples
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		def samples = Samples.fromParams(params,{sheet,schema -> samplesheetToList(sheet, schema)})
		def sr_ch = samples
				.filter({it.reads_short_1})
				.map({[it.subMap('sample_id'),[it.reads_short_1,it.short_reads_2].findAll({it})]})
		def lr_ch = samples
			.filter({it.reads_long})
			.map({[it.subMap('sample_id'),it.reads_long]})

		// Filter reads
		READS_QC(sr_ch,lr_ch)
		sr_ch = READS_QC.out.short_filtered
		lr_ch = READS_QC.out.long_filtered
		
		ASSEMBLE(params.assembler.name,sr_ch,lr_ch)
		ASSEMBLIES_QC(ASSEMBLE.out.fasta,sr_ch,lr_ch)
		AMR_ANNOT(params.amr_annot,ASSEMBLE.out.fasta,sr_ch,lr_ch)
}

