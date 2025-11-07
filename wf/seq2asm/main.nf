
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

		// Short reads only assemblies
		SHORT_SPADES(fqs_ch.filter({opts.short_spades}),Channel.empty())
		SHORT_UNICYCLER(fqs_ch.filter({opts.short_unicycler}),Channel.empty())
		
		// Long reads only assemblies
		LONG_HYBRACTER(Channel.empty(),fql_ch.filter({opts.long_hybracter}))
		LONG_UNICYCLER(Channel.empty(),fql_ch.filter({opts.long_unicycler}))
		LONG_FLYE_MEDAKA(Channel.empty(),fql_ch.filter({opts.long_flye_medaka}))

		// Hybrid assemblies
		HYBRID_HYBRACTER(fqs_ch.filter({opts.hybrid_hybracter}),fql_ch.filter({opts.hybrid_hybracter}))
		HYBRID_UNICYCLER(fqs_ch.filter({opts.hybrid_unicycler}),fql_ch.filter({opts.hybrid_unicycler}))
		HYBRID_FLYE_MEDAKA_PILON(fqs_ch.filter({opts.hybrid_flye_medaka_pilon}),fql_ch.filter({opts.hybrid_flye_medaka_pilon}))
		
	emit:
		fasta = Channel.empty().mix(
			SHORT_SPADES.out.fasta.map({meta,x -> [meta+[assembly_name:'short_spades'],x]}),
			SHORT_UNICYCLER.out.fasta.map({meta,x -> [meta+[assembly_name:'short_unicycler'],x]}),
			LONG_FLYE_MEDAKA.out.fasta.map({meta,x -> [meta+[assembly_name:'long_flye_medaka'],x]}),
			LONG_UNICYCLER.out.fasta.map({meta,x -> [meta+[assembly_name:'long_unicycler'],x]}),
			LONG_HYBRACTER.out.fasta.map({meta,x -> [meta+[assembly_name:'long_hybracter'],x]}),
		  HYBRID_UNICYCLER.out.fasta.map({meta,x -> [meta+[assembly_name:'hybrid_unicycler'],x]}),
		  HYBRID_HYBRACTER.out.fasta.map({meta,x -> [meta+[assembly_name:'hybrid_hybracter'],x]}),
		  HYBRID_FLYE_MEDAKA_PILON.out.fasta.map({meta,x -> [meta+[assembly_name:'hybrid_flye_medaka_pilon'],x]})
		)
		dir = Channel.empty().mix(
			SHORT_SPADES.out.dir.map({meta,x -> [meta+[assembly_name:'short_spades'],x]}),
			SHORT_UNICYCLER.out.dir.map({meta,x -> [meta+[assembly_name:'short_unicycler'],x]}),
			LONG_FLYE_MEDAKA.out.dir.map({meta,x -> [meta+[assembly_name:'long_flye_medaka'],x]}),
			LONG_UNICYCLER.out.dir.map({meta,x -> [meta+[assembly_name:'long_unicycler'],x]}),
			LONG_HYBRACTER.out.dir.map({meta,x -> [meta+[assembly_name:'long_hybracter'],x]}),
		  HYBRID_UNICYCLER.out.dir.map({meta,x -> [meta+[assembly_name:'hybrid_unicycler'],x]}),
		  HYBRID_HYBRACTER.out.dir.map({meta,x -> [meta+[assembly_name:'hybrid_hybracter'],x]}),
		  HYBRID_FLYE_MEDAKA_PILON.out.dir.map({meta,x -> [meta+[assembly_name:'hybrid_flye_medaka_pilon'],x]})
		)
}




// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.samplesheet = null
params.long_reads  = []
params.short_reads = []
params.long_unicycler           = false
params.long_hybracter           = false
params.long_flye_medaka         = false
params.short_spades             = false
params.short_unicycler          = false
params.hybrid_unicycler         = false
params.hybrid_hybracter         = false
params.hybrid_flye_medaka_pilon = false


import AmrUtils


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Extract Long_read and Short_read channels from params 
		def ss = AmrUtils.parse_generic_params(params,samplesheetToList)

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		ss.fqs_ch = ss.fqs_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		ss.fqs_ch = ss.fqs_ch.fq.mix(SAMTOOLS_FASTQ(ss.fqs_ch.bam))

		// Run assemblers
		SEQ2ASM(params,ss.fqs_ch,ss.fql_ch)

	publish:
		fasta            = SEQ2ASM.out.fasta
		assembler_output = SEQ2ASM.out.dir
}

output {
	fasta {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/assembly.fasta"}
	}
	assembler_output {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/assembler"}
	}
}













