
include { ORGFINDER_DB_DOWNLOAD } from '../subworkflows/orgfinder/db/main.nf'

workflow {
	main:
		def db = ORGFINDER_DB_DOWNLOAD()
		
}

