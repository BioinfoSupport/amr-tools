

include {READS_QC}    from './wf/reads_qc/main.nf'
include {ASSEMBLE}   from './wf/assemble/main.nf'
//include {ASM_QC}    from './asm_qc/main.nf'
//include {AMR_ANNOT} from './amr_annot/main.nf'



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
			Samples.short_reads_channel(samples),
			Samples.long_reads_channel(samples)
		)

/*		
		samples
		.filter({it.reads_long})
		.map({[it.subMab('sample_id','assembler_name','assembler_args'),it.reads_long]})
*/


		if (params.mode=="assembly_qc") {
			log.info("assemble not implemented yet !")
		}
		if (params.mode=="amr_annot") {
			log.info("assemble not implemented yet !")
		}
		
}

