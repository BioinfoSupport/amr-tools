
include { ORGFINDER_DB } from './subworkflows/orgfinder_db/main.nf'
include { FASTANI } from './modules/fastani/main.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		ORGFINDER_DB()

		def query_ch = Channel.fromPath(params.query).map({
				def id = it.name.replaceAll(/\.(fasta|fna|fa)$/,'')
				[[sample_id:id],it]
		})

		def ref_ch = Channel.fromPath(params.ref).collect().map({[it]})
		query_ch
				.combine(ref_ch)
				.map({m, q, r -> tuple(m, q, r)})
				| FASTANI

	publish:
		fastani_tsv = FASTANI.out
		db = ORGFINDER_DB.out
}

output {
	db {
		path "db/"
	}
	fastani_tsv {
		path { m,x -> x >> "${m.sample_id}.fastani"}
	}
}




