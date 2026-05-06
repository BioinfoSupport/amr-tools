


include { PLING_CLUSTER_ALIGN   } from '../../modules/pling/cluster_align.nf'
workflow TEST_PLING {
	main:
		Channel.fromPath("${moduleDir}/../../data/plasmids/fastas/*.fna")
		| collectFile(name: 'plasmids.fasta')
		| map {["all_plasmids",it]}
		| PLING_CLUSTER_ALIGN
}


include { ORGFINDER_DB_DOWNLOAD } from '../../subworkflows/orgfinder/db/main.nf'
include { ORGFINDER_DETECT      } from '../../subworkflows/orgfinder/detect/main.nf'
workflow TEST_ORGFINDER {
	main:
		def fa_ch = Channel.fromPath("${moduleDir}/assets/*.fasta")
		ORGFINDER_DB_DOWNLOAD()
		ORGFINDER_DETECT(fa_ch,ORGFINDER_DB_DOWNLOAD.out)
}


workflow {
	//TEST_PLING()
	TEST_ORGFINDER()
}