

include {READS_QC}    from './wf/reads_qc/main.nf'
//include {ASSEMBLE}   from './wf/assemble/main.nf'
//include {ASM_QC}    from './asm_qc/main.nf'
//include {AMR_ANNOT} from './amr_annot/main.nf'



// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
import Samples
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.mode = null
params.csv = null
params.reads = [long:null,short:[]]
params.assembler = [name:null,args:null]
params.assemblies = [fasta:null,info:null]
params.reads_qc = [
		limit_long_reads_len: 1000000000,
		limit_short_reads_len: 500000000
]

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		def samples = Samples.fromParams(params,{sheet,schema -> samplesheetToList(sheet, schema)})
		
		if (params.mode=="reads_qc") {
			READS_QC(params.reads_qc,Samples.short_reads_channel(samples),Samples.long_reads_channel(samples))	
		} else if (params.mode=="assemble") {
			log.info("assemble not implemented yet !")
		}
		
}

