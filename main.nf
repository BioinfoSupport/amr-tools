

include {READS_QC}  from './wf/reads_qc/main.nf'
include {ASSEMBLE}  from './wf/assemble/main.nf'
//include {ASM_QC}    from './wf/asm_qc/main.nf'
//include {AMR_ANNOT} from './wf/amr_annot/main.nf'



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
		READS_QC(
			Samples.short_reads_channel(samples),
			Samples.long_reads_channel(samples)
		)
		ASSEMBLE(
			params.assembler.name,
			READS_QC.out.short_filtered,
			READS_QC.out.long_filtered
		)
		
		
}

