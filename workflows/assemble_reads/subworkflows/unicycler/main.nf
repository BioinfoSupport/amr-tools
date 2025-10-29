
include { UNICYCLER as UNICYCLER_RUN  } from './modules/unicycler'
include { UNICYCLER_ADAPT             } from './modules/unicycler/adapt'

workflow UNICYCLER {
	take:
		fqs_ch
	  fql_ch
	main:
		UNICYCLER_RUN(fqs_ch.join(fql_ch,remainder:true).map({meta,fqs,fql -> [meta,fqs?:[],fql?:[]]})) | UNICYCLER_ADAPT
	emit:
		fasta = UNICYCLER_ADAPT.out.fasta
		dir   = UNICYCLER_RUN.out
}
