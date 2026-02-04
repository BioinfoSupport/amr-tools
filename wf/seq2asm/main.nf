
include { SAMTOOLS_FASTQ                                 } from './modules/samtools/fastq'
include { UNICYCLER as LONG_UNICYCLER                    } from './subworkflows/unicycler'
include { UNICYCLER as SHORT_UNICYCLER                   } from './subworkflows/unicycler'
include { UNICYCLER as HYBRID_UNICYCLER                  } from './subworkflows/unicycler'
include { HYBRACTER as LONG_HYBRACTER                    } from './subworkflows/hybracter'
include { HYBRACTER as HYBRID_HYBRACTER                  } from './subworkflows/hybracter'
include { SPADES    as SHORT_SPADES                      } from './subworkflows/spades'
include { FLYE_MEDAKA_PILON as HYBRID_FLYE_MEDAKA_PILON  } from './subworkflows/flye_medaka_pilon'
include { FLYE_MEDAKA_PILON as LONG_FLYE_MEDAKA          } from './subworkflows/flye_medaka_pilon'



workflow SEQ2ASM {
	take:
		opts
		fqs_ch    // channel: [ val(meta), path(short_reads) ]
		fql_ch    // channel: [ val(meta), path(long_reads) ]
	main:

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		fql_ch = fql_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		fql_ch = fql_ch.fq.mix(SAMTOOLS_FASTQ(fql_ch.bam))

		
		def assemblies_fasta = Channel.empty()
		def assemblies_dir = Channel.empty()
		
		// Short reads only assemblies
		if (opts.name=="short_spades") {
			SHORT_SPADES(fqs_ch,Channel.empty())
			assemblies_fasta = assemblies_fasta.mix(SHORT_SPADES.out.fasta.map({meta,x -> [meta+[assembly_name:'short_spades'],x]}))
			assemblies_dir = assemblies_dir.mix(SHORT_SPADES.out.dir.map({meta,x -> [meta+[assembly_name:'short_spades'],x]}))
		}
		if (opts.name=="short_unicycler") {
			SHORT_UNICYCLER(fqs_ch,Channel.empty())
			assemblies_fasta = assemblies_fasta.mix(SHORT_UNICYCLER.out.fasta.map({meta,x -> [meta+[assembly_name:'short_unicycler'],x]}))
			assemblies_dir = assemblies_dir.mix(SHORT_UNICYCLER.out.dir.map({meta,x -> [meta+[assembly_name:'short_unicycler'],x]}))
		}
		
		// Long reads only assemblies
		if (opts.name=="long_hybracter") {
			LONG_HYBRACTER(Channel.empty(),fql_ch)
			assemblies_fasta = assemblies_fasta.mix(LONG_HYBRACTER.out.fasta.map({meta,x -> [meta+[assembly_name:'long_hybracter'],x]}))
			assemblies_dir = assemblies_dir.mix(LONG_HYBRACTER.out.dir.map({meta,x -> [meta+[assembly_name:'long_hybracter'],x]}))
		}
		if (opts.name=="long_unicycler") {
			LONG_UNICYCLER(Channel.empty(),fql_ch)
			assemblies_fasta = assemblies_fasta.mix(LONG_UNICYCLER.out.fasta.map({meta,x -> [meta+[assembly_name:'long_unicycler'],x]}))
			assemblies_dir = assemblies_dir.mix(LONG_UNICYCLER.out.dir.map({meta,x -> [meta+[assembly_name:'long_unicycler'],x]}))
		}
		if (opts.name=="long_flye_medaka") {
			LONG_FLYE_MEDAKA(Channel.empty(),fql_ch)
			assemblies_fasta = assemblies_fasta.mix(LONG_FLYE_MEDAKA.out.fasta.map({meta,x -> [meta+[assembly_name:'long_flye_medaka'],x]}))
			assemblies_dir = assemblies_dir.mix(LONG_FLYE_MEDAKA.out.dir.map({meta,x -> [meta+[assembly_name:'long_flye_medaka'],x]}))
		}

		// Hybrid assemblies
		if (opts.name=="hybrid_hybracter") {
			HYBRID_HYBRACTER(fqs_ch,fql_ch)
			assemblies_fasta = assemblies_fasta.mix(HYBRID_HYBRACTER.out.fasta.map({meta,x -> [meta+[assembly_name:'hybrid_hybracter'],x]}))
			assemblies_dir = assemblies_dir.mix(HYBRID_HYBRACTER.out.dir.map({meta,x -> [meta+[assembly_name:'hybrid_hybracter'],x]}))
		}
		if (opts.name=="hybrid_unicycler") {
			HYBRID_UNICYCLER(fqs_ch,fql_ch)
			assemblies_fasta = assemblies_fasta.mix(HYBRID_UNICYCLER.out.fasta.map({meta,x -> [meta+[assembly_name:'hybrid_unicycler'],x]}))
			assemblies_dir = assemblies_dir.mix(HYBRID_UNICYCLER.out.dir.map({meta,x -> [meta+[assembly_name:'hybrid_unicycler'],x]}))
		}
		if (opts.name=="hybrid_flye_medaka_pilon") {
			HYBRID_FLYE_MEDAKA_PILON(fqs_ch,fql_ch)
			assemblies_fasta = assemblies_fasta.mix(HYBRID_FLYE_MEDAKA_PILON.out.fasta.map({meta,x -> [meta+[assembly_name:'hybrid_flye_medaka_pilon'],x]}))
			assemblies_dir = assemblies_dir.mix(HYBRID_FLYE_MEDAKA_PILON.out.dir.map({meta,x -> [meta+[assembly_name:'hybrid_flye_medaka_pilon'],x]}))
		}
		
	emit:
		fasta = assemblies_fasta
		dir = assemblies_dir
}




// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
import Readsets
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.readsets = [
	csv         : null,
	long_reads  : [],
	short_reads : []
]

params.assembler = [
	name: "long_flye_medaka",
	args: ""
]

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Extract Long_read and Short_read channels from params
		def readsets = Readsets.fromParams(params.readsets,{sheet,schema -> samplesheetToList(sheet, schema)})
		
		// Run assemblers
		SEQ2ASM(params.assembler,readsets.short_reads,readsets.long_reads)

		assemblies = SEQ2ASM.out.fasta
			.join(SEQ2ASM.out.dir,remainder:true)
			.map({m,fa,dir -> m + [assembly_fasta:fa,assembler_output:dir]})

	publish:
		assemblies   = assemblies
		readsets_csv = readsets.flat_csv()
}

output {
	readsets_csv {
    index {
    	path 'indexes/seq2asm_input_readsets.csv'
    	header true
    }
	}
	assemblies {
    path { x ->
    	x.assembly_fasta >> "samples/${x.sample_id}/assemblies/${x.assembly_name}/asm.fasta"
    	x.assembler_output >> "samples/${x.sample_id}/assemblies/${x.assembly_name}/assembler_output"
    }
    index {
    	path 'indexes/seq2asm_assemblies.csv'
    	header true
    }
	}
}













