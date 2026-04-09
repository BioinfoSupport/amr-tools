

include { FASTANI } from './modules/fastani/main.nf'


process ORGFINDER_DB_ADAPT {
	input:
		path('db')
	output:
		path('db/fna/*.fna'), emit: fna
	script:
	"""
	"""
}


workflow ORGFINDER_DETECT {
	take:
		fa_ch   // channel: tuple val(meta), path(fasta)
		db_ch   // channel: path(db)
	main:
		def query_ch = fa_ch.collect().map({[it]})
		def ref_ch = ORGFINDER_DB_ADAPT(db_ch).fna.collect().map({[it]})
		FASTANI(
			query_ch
					.combine(ref_ch)
					.map({q, r -> tuple("orgfinder", q, r)})
		)
	emit:
		res = FASTANI.out
}
