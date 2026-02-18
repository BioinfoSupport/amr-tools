
include { SAMTOOLS_FASTQ as CONVERT_LONG_BAM_TO_FASTQ } from './modules/samtools/fastq' 
include { READS_FILTER                                } from './wf/reads_filter/main.nf'
include { READS_QC as FILTERED_READS_QC               } from './wf/reads_qc/main.nf'
include { ASSEMBLE                                    } from './wf/assemble/main.nf'
include { ASSEMBLIES_QC                               } from './wf/assemblies_qc/main.nf'
include { AMR_ANNOT                                   } from './wf/amr_annot/main.nf'
include { IDENTITY as COPY_ASSEMBLY_FASTA             } from './modules/identity/main.nf'
include { IDENTITY as COPY_SHORT_READS                } from './modules/identity/main.nf'
include { IDENTITY as COPY_LONG_READS                 } from './modules/identity/main.nf'

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
				.map({[it.subMap('sample_id'),[it.reads_short_1,it.reads_short_2].findAll({it})]})
				| COPY_SHORT_READS
		
		def lr_ch = samples
			.filter({it.reads_long})
			.map({[it.subMap('sample_id'),it.reads_long]})
			| COPY_LONG_READS

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		lr_ch = lr_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		lr_ch = lr_ch.fq.mix(CONVERT_LONG_BAM_TO_FASTQ(lr_ch.bam))

		// Filter reads
		READS_FILTER(sr_ch,lr_ch)
		FILTERED_READS_QC(READS_FILTER.out.short_filtered,READS_FILTER.out.long_filtered)
		
		// Run denovo assembly process when 
		// Determine assemblies that are done and thus that have to be run
		def asm_ch = samples.branch({
			done: it.assembly_fasta
			todo: true
		})
		asm_ch.done.map({[it.subMap('sample_id'),it.assembly_fasta]}) 
		| COPY_ASSEMBLY_FASTA
		ASSEMBLE(
				params.assembler.name,
				READS_FILTER.out.short_filtered.join(asm_ch.todo.map({[it.subMap('sample_id')]})),
				READS_FILTER.out.long_filtered.join(asm_ch.todo.map({[it.subMap('sample_id')]}))
		)
		def asm_fa_ch = Channel.empty().mix(
				COPY_ASSEMBLY_FASTA.out.map({[[it[0],[assembler_name:'none']],it[1]]}),
				ASSEMBLE.out.fasta.map({[[it[0],[assembler_name:params.assembler.name]],it[1]]})
		)

		// Compute assemblies QC stats
		ASSEMBLIES_QC(
			asm_fa_ch,
			READS_FILTER.out.short_filtered.combine(asm_fa_ch.map({m,x->m}),by:0).map({[[it[0],it[2]],it[1]]}),
			READS_FILTER.out.long_filtered.combine(asm_fa_ch.map({m,x->m}),by:0).map({[[it[0],it[2]],it[1]]})
		)

		

		// Finally run AMR annotations
		AMR_ANNOT(
			params.amr,
			asm_fa_ch,
			READS_FILTER.out.short_filtered.combine(asm_fa_ch.map({m,x->m}),by:0).map({[[it[0],it[2]],it[1]]}),
			READS_FILTER.out.long_filtered.combine(asm_fa_ch.map({m,x->m}),by:0).map({[[it[0],it[2]],it[1]]})
		)

	publish:
		reads_long_filtered = READS_FILTER.out.long_filtered
		reads_short_filtered = READS_FILTER.out.short_filtered
		filtered_reads_multiqc_html = FILTERED_READS_QC.out.multiqc_html
		assemblies_fasta = asm_fa_ch
		assemblies_multiqc_txt = ASSEMBLIES_QC.out.assembly_multiqc_txt
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
	filtered_reads_multiqc_html {
		path { x -> x >> "filtered_reads_multiqc.html"}
	}
	
	assemblies_fasta {
		path { m,x -> x >> "assemblies/${m[1].assembler_name}/${m[0].sample_id}.fasta"}
	}
	assemblies_multiqc_txt {
		path { x -> x >> "assemblies_multiqc.txt"}
	}
	
/*
	orgfinder {
		path { m,x -> x >> "assemblies/${m.assembler_name}/${m.sample_id}/orgfinder"}
	}
*/
}







