

include {READS_QC}    from './wf/reads_qc/main.nf'
//include {ASSEMBLE}   from './wf/assemble/main.nf'
//include {ASM_QC}    from './asm_qc/main.nf'
//include {AMR_ANNOT} from './amr_annot/main.nf'



// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
import Samples
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params = [
	mode      : null,
	csv       : null,
	reads     : [long:null,short:[]],
	assembler : [name:null,args:null],
	assemblies: [fasta:null,info:null],
	reads_qc: [
		limit_long_reads_len: 1000000000,
		limit_short_reads_len: 1000000000
	]
]

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		def samples = Samples.fromParams(params,{sheet,schema -> samplesheetToList(sheet, schema)})
		
		if (params.mode=="reads_qc") {
			READS_QC(params,Readsets.short_reads_channel(readsets),Readsets.long_reads_channel(readsets))	
		} else if (params.mode=="assemble") {
			log.info("assemble not implemented yet !")
		}
		
}

