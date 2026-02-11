

include { SAMTOOLS_FASTQ                           } from './modules/samtools/fastq'
include { FLYE_MEDAKA as LONG_FLYE_MEDAKA          } from './subworkflows/flye_medaka_pilon'

/*
include { FLYE_MEDAKA_PILON as HYBRID_FLYE_MEDAKA_PILON  } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER as LONG_UNICYCLER                    } from './subworkflows/unicycler'
include { UNICYCLER as SHORT_UNICYCLER                   } from './subworkflows/unicycler'
include { UNICYCLER as HYBRID_UNICYCLER                  } from './subworkflows/unicycler'
include { HYBRACTER as LONG_HYBRACTER                    } from './subworkflows/hybracter'
include { HYBRACTER as HYBRID_HYBRACTER                  } from './subworkflows/hybracter'
include { SPADES    as SHORT_SPADES                      } from './subworkflows/spades'
*/




workflow ASSEMBLE {
	take:
		assembler_name
		fqs_ch    // channel: [ val(meta), path(short_reads) ]
		fql_ch    // channel: [ val(meta), path(long_reads) ]
	main:

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format if needed
		fql_ch = fql_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		fql_ch = fql_ch.fq.mix(SAMTOOLS_FASTQ(fql_ch.bam))
		
		def assemblies = [
			fasta: Channel.empty(),
			dir: Channel.empty()
		]

		switch(assembler_name) {
			case 'long_flye_medaka':
				LONG_FLYE_MEDAKA(Channel.empty(),fql_ch)
				assemblies.fasta = assemblies.fasta.mix(LONG_FLYE_MEDAKA.out.fasta.map({meta,x -> [meta+[assembler_name:'long_flye_medaka'],x]}))
				assemblies.dir = assemblies.dir.mix(LONG_FLYE_MEDAKA.out.dir.map({meta,x -> [meta+[assembler_name:'long_flye_medaka'],x]}))
				break
			default:
				error "Unknown assembler name :${asssembler_name}"
		}

	emit:
		fasta = assemblies.fasta
		dir = assemblies.dir
}

