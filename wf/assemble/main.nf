

include { SAMTOOLS_FASTQ                } from './modules/samtools/fastq'
include { SAMTOOLS_FAIDX                } from './modules/samtools/faidx'
include { LONG_FLYE                     } from './subworkflows/flye_medaka_pilon'
include { LONG_FLYE_MEDAKA              } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER as LONG_UNICYCLER   } from './subworkflows/unicycler'
include { HYBRACTER as LONG_HYBRACTER   } from './subworkflows/hybracter'
include { UNICYCLER as SHORT_UNICYCLER  } from './subworkflows/unicycler'
include { SPADES    as SHORT_SPADES     } from './subworkflows/spades'
include { UNICYCLER as HYBRID_UNICYCLER } from './subworkflows/unicycler'
include { HYBRACTER as HYBRID_HYBRACTER } from './subworkflows/hybracter'
include { HYBRID_FLYE_MEDAKA_PILON      } from './subworkflows/flye_medaka_pilon'


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
		
		def assemblies = [fasta:Channel.empty(),dir:Channel.empty()]
		switch(assembler_name) {
			
			// Long reads assemblers
			case 'long_flye':
				LONG_FLYE(fql_ch)
				assemblies = [fasta:LONG_FLYE.out.fasta,dir:LONG_FLYE.out.dir]
				break
			case 'long_flye_medaka':
				LONG_FLYE_MEDAKA(fql_ch)
				assemblies = [fasta:LONG_FLYE_MEDAKA.out.fasta,dir:LONG_FLYE_MEDAKA.out.dir]
				break
			case 'long_unicycler':
				LONG_UNICYCLER(Channel.empty(),fql_ch)
				assemblies = [fasta:LONG_UNICYCLER.out.fasta,dir:LONG_UNICYCLER.out.dir]
				break
			case 'long_hybracter':
				LONG_HYBRACTER(Channel.empty(),fql_ch)
				assemblies = [fasta:LONG_HYBRACTER.out.fasta,dir:LONG_HYBRACTER.out.dir]
				break
				
			// Short reads assemblers
			case 'short_unicycler':
				SHORT_UNICYCLER(fqs_ch,Channel.empty())
				assemblies = [fasta:SHORT_UNICYCLER.out.fasta,dir:SHORT_UNICYCLER.out.dir]
				break
			case 'short_spades':
				SHORT_SPADES(fqs_ch,Channel.empty())
				assemblies = [fasta:SHORT_SPADES.out.fasta,dir:SHORT_SPADES.out.dir]
				break
				
			// Hybrid assemblers
			case 'hybrid_unicycler':
				HYBRID_UNICYCLER(fqs_ch,fql_ch)
				assemblies = [fasta:HYBRID_UNICYCLER.out.fasta,dir:HYBRID_UNICYCLER.out.dir]
				break
			case 'hybrid_hybracter':
				HYBRID_HYBRACTER(fqs_ch,fql_ch)
				assemblies = [fasta:HYBRID_HYBRACTER.out.fasta,dir:HYBRID_HYBRACTER.out.dir]
				break
			case 'hybrid_flye_medaka_pilon':
				HYBRID_FLYE_MEDAKA_PILON(fqs_ch,fql_ch)
				assemblies = [fasta:HYBRID_FLYE_MEDAKA_PILON.out.fasta,dir:HYBRID_FLYE_MEDAKA_PILON.out.dir]
				break
			default:
				error "Unknown assembler name :${assembler_name}"
		}
		
		SAMTOOLS_FAIDX(assemblies.fasta)
		
	emit:
		fasta = assemblies.fasta
		fai   = SAMTOOLS_FAIDX.out
		dir = assemblies.dir
}

