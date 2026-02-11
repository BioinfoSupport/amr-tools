
include { SPADES as SPADES_RUN } from './modules/spades'
include { SPADES_ADAPT         } from './modules/spades/adapt'

workflow SPADES {
	take:
		fqs_ch
		fql_ch
	main:
		SPADES_RUN(
			fqs_ch
			.join(fql_ch,remainder:true)
			.map({meta,fqs,fql -> [meta,fqs?:[],fql?:[]]})
		) 
		| SPADES_ADAPT
	emit:
		fasta = SPADES_ADAPT.out.fasta
		dir   = SPADES_RUN.out
}
