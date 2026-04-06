#!/usr/bin/env nextflow

include { ORGFINDER_DB_DOWNLOAD } from './subworkflows/orgfinder/mains.nf'

workflow TESTS {
	main:
		ORGFINDER_DB_DOWNLOAD()
}














