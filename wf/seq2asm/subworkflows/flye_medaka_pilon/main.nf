
include { FLYE             } from './modules/flye'
include { FLYE_ADAPT       } from './modules/flye/adapt'
include { MEDAKA_CONSENSUS } from './modules/medaka/consensus'
include { MEDAKA_ADAPT     } from './modules/medaka/adapt'
include { ORGANIZE_FILES   } from './modules/organize_files'
include { PILON_POLISH as PILON_POLISH_ROUND1 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND2 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND3 } from './subworkflows/pilon_polish'


workflow FLYE_MEDAKA_PILON {
	take:
		fqs_ch
		fql_ch
	main:
		FLYE(fql_ch) | FLYE_ADAPT
		FLYE_ADAPT.out.fasta.join(fql_ch) | MEDAKA_CONSENSUS | MEDAKA_ADAPT
		PILON_POLISH_ROUND1(MEDAKA_ADAPT.out.fasta,fqs_ch)
		PILON_POLISH_ROUND2(PILON_POLISH_ROUND1.out.fasta,fqs_ch)
		PILON_POLISH_ROUND3(PILON_POLISH_ROUND2.out.fasta,fqs_ch)
		ORGANIZE_FILES(
			FLYE.out
			.join(MEDAKA_CONSENSUS.out, remainder: true)
			.join(PILON_POLISH_ROUND1.out.dir, remainder: true)
			.join(PILON_POLISH_ROUND2.out.dir, remainder: true)
			.join(PILON_POLISH_ROUND3.out.dir, remainder: true)
			.map({meta,x1,x2,x3,x4,x5 -> tuple(meta,[[x1,'01_flye'],[x2,'02_medaka'],[x3,'03_pilon_round1'],[x4,'04_pilon_round2'],[x5,'05_pilon_round3']].findAll({x,y -> x}))})
		)
		
		fa_ch = MEDAKA_ADAPT.out.fasta
		.join(PILON_POLISH_ROUND3.out.fasta,remainder: true)
		.map({meta,fa1,fa2 -> [meta,fa2?fa2:fa1]})
	emit:
		fasta = fa_ch
		dir = ORGANIZE_FILES.out
}
