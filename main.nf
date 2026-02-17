

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


		def long_filtered_ch = Channel.empty()
		def short_filtered_ch = Channel.empty()
		def assemblies_fasta_ch = Channel.empty()
		def assemblies_fai_ch = Channel.empty()
		def orgfinder_ch = Channel.empty()
		switch(params.subcommand) {
			case 'reads_filter':
				READS_FILTER(sr_ch,lr_ch)
				long_filtered_ch = READS_FILTER.out.long_filtered
				short_filtered_ch = READS_FILTER.out.short_filtered
				break
			case 'assemble':
				ASSEMBLE(params.assembler.name,sr_ch,lr_ch)
				break
			case 'assembly_qc':
				ASSEMBLIES_QC(ASSEMBLE.out.fasta,sr_ch,lr_ch)
				assemblies_fasta_ch = ASSEMBLE.out.fasta
				assemblies_fai_ch = ASSEMBLE.out.fasta
				break
			case 'amr_annot':
				AMR_ANNOT(params.amr_annot,ASSEMBLE.out.fasta,sr_ch,lr_ch)
				orgfinder_ch = AMR_ANNOT.out.orgfinder
				break
			default:
				READS_FILTER(sr_ch,lr_ch)
				ASSEMBLE(
					params.assembler.name,
					READS_FILTER.out.short_filtered,
					READS_FILTER.out.long_filtered
				)
				AMR_ANNOT(params.amr_annot,ASSEMBLE.out.fasta,READS_FILTER.out.short_filtered,READS_FILTER.out.long_filtered)
				long_filtered_ch = READS_FILTER.out.long_filtered
				short_filtered_ch = READS_FILTER.out.short_filtered
				assemblies_fasta_ch = ASSEMBLE.out.fasta
				assemblies_fai_ch = ASSEMBLE.out.fasta
				orgfinder_ch = AMR_ANNOT.out.orgfinder
		}		
		
	publish:
		reads_long_filtered = long_filtered_ch
		reads_short_filtered = short_filtered_ch
		assemblies_fasta = assemblies_fasta_ch
		assemblies_fai = assemblies_fai_ch
		orgfinder = orgfinder_ch
}

output {
	
	// -------------------
  // reads_filter command
  // -------------------
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



  // -------------------
  // assemble command
  // -------------------
	assemblies_fasta {
		path { m,x -> x >> "assemblies/${params.assembler.name}/${m.sample_id}.fasta"}
	}
	assemblies_fai {
		path { m,x -> x >> "assemblies/${params.assembler.name}/${m.sample_id}.fasta.fai"}
	}
	
	
	
	
  // -------------------
  // amr_annot command
  // -------------------	
	orgfinder {
		path { m,x -> x >> "assemblies/${params.assembler.name}/${m.sample_id}/orgfinder"}
	}
}







