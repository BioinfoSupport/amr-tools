

include {READS_FILTER}  from './wf/reads_filter/main.nf'
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
		switch(params.command) {
			case 'reads_filter':
				READS_FILTER(sr_ch,lr_ch)
				break
			case 'assemble':
				ASSEMBLE(params.assembler.name,sr_ch,lr_ch)
				break
			case 'assembly_qc':
				ASSEMBLIES_QC(ASSEMBLE.out.fasta,sr_ch,lr_ch)
				break
			case 'amr_annot':
				AMR_ANNOT(params.amr_annot,ASSEMBLE.out.fasta,sr_ch,lr_ch)
				break
			default:
				error "Unknown command :${params.command}"
		}		
		
	publish:
		reads_long_filtered = (params.command in ['reads_filter'])?READS_FILTER.out.long_filtered:Channel.empty()
		reads_short_filtered = (params.command in ['reads_filter'])?READS_FILTER.out.short_filtered:Channel.empty()
		assemblies_fasta = (params.command in ['assemble'])?ASSEMBLE.out.fasta:Channel.empty()
		orgfinder = (params.command in ['amr_annot'])?AMR_ANNOT.out.orgfinder:Channel.empty()
}

output {
	reads_long_filtered {
		path { m,x -> x >> "reads/filtered_long/${m.sample_id}.fastq.gz"}
	}
	
	reads_short_filtered {
		path { m,x -> 
			x = x instanceof List?x:[x]
			x.eachWithIndex { xi, i ->
      	xi >> "reads/filtered_short/${m.sample_id}_R${i+1}.fastq.gz"
      }
		}
	}

	assemblies_fasta {
		path { m,x -> x >> "assemblies/${params.assembler.name}/${m.sample_id}/asm.fasta"}
	}

	orgfinder {
		path { m,x -> x >> "assemblies/${params.assembler.name}/${m.sample_id}/orgfinder"}
	}
}







