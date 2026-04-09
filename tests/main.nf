
include { ORGFINDER_DB_DOWNLOAD } from '../subworkflows/orgfinder/db/main.nf'
include { ORGFINDER_DETECT      } from '../subworkflows/orgfinder/detect/main.nf'

workflow {
	main:
		def fa_ch = Channel.fromPath("${moduleDir}/assets/*.fasta")
		ORGFINDER_DB_DOWNLOAD()
		ORGFINDER_DETECT(fa_ch,ORGFINDER_DB_DOWNLOAD.out)
}

