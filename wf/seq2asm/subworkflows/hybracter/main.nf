
include { HYBRACTER as HYBRACTER_RUN  } from './modules/hybracter'
include { HYBRACTER_ADAPT             } from './modules/hybracter/adapt'

workflow HYBRACTER {
	take:
		fqs_ch
	  fql_ch
	  args
	main:
		HYBRACTER_RUN(
			fqs_ch
			.join(fql_ch,remainder:true)
			.map({meta,fqs,fql -> [meta,fqs?:[],fql?:[]]})
			.combine(Channel.from(args?:""))
		) 
		| HYBRACTER_ADAPT
	emit:
		fasta = HYBRACTER_ADAPT.out.fasta
		dir   = HYBRACTER_RUN.out
}
