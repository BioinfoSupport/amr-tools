

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

		// Parse parameters and get a channel of samples objects
		def samples = Samples.fromParams(params,{sheet,schema -> samplesheetToList(sheet, schema)})
		
		// Extract reads channels from samples channel
		def sr_ch = samples
				.filter({it.reads_short_1})
				.map({[it.subMap('sample_id'),[it.reads_short_1,it.short_reads_2].findAll({it})]})
		def lr_ch = samples
			.filter({it.reads_long})
			.map({[it.subMap('sample_id'),it.reads_long]})

		// Filter reads
		READS_FILTER(sr_ch,lr_ch)
		
		// Determine assemblies that are done and thus that have to be run
		def asm_ch = samples.branch({
			done: it.assembly_fasta
			todo: true
		})
		ASSEMBLE(
				params.assembler.name,
				sr_ch.join(asm_ch.todo.map({[it.subMap('sample_id')]})),
				lr_ch.join(asm_ch.todo.map({[it.subMap('sample_id')]}))
		)
		def asm_fa_ch = Channel.empty().mix(
			//asm_ch.done.map({[it.subMap('sample_id'),[assembler_name:'none'], file(it.assembly_fasta)]}),
			ASSEMBLE.out.fasta.map({[it[0],[assembler_name:params.assembler.name],it[1]]})
		)

		// Finally run AMR annotations
		AMR_ANNOT(
			params.amr_annot,
			asm_fa_ch.map({[it[0]+it[1],it[2]]}),
			asm_fa_ch.join(READS_FILTER.out.short_filtered).map({[it[0]+it[1]?:[:],it[3]]}),
			asm_fa_ch.join(READS_FILTER.out.long_filtered).map({[it[0]+it[1]?:[:],it[3]]})
		)

//asm_fa_ch.map({ m,m2,x -> "assemblies/${m2.assembler_name}/${m.sample_id}.fasta"}).view()

	publish:
		reads_long_filtered = READS_FILTER.out.long_filtered
		reads_short_filtered = READS_FILTER.out.short_filtered
		assemblies_fasta = asm_fa_ch
		//orgfinder = AMR_ANNOT.out.orgfinder
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
		path { m,m2,x -> x >> "assemblies/${m2.assembler_name}/${m.sample_id}.fasta"}
	}
/*
	orgfinder {
		path { m,x -> x >> "assemblies/${m.assembler_name}/${m.sample_id}/orgfinder"}
	}
*/
}







