

include { NCBI_DATASET_DOWNLOAD_GENOME } from './modules/ncbi/dataset/main.nf'
include { NCBI_TAXDUMP_DOWNLOAD        } from './modules/ncbi/taxdump/main.nf'
include { RSCRIPT                      } from './modules/rscript/main.nf'


process GENOMES_AGGREGATE {
    container 'docker.io/staphb/ncbi-datasets:18.18.0'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
  		path('genomes/genome*')
  	output:
  		path('genomes/')
    script:
	    """
			# Make the tsv file with all accession numbers
			cat genomes/*/ncbi_dataset/data/assembly_data_report.jsonl \
			  | dataformat tsv genome --force --fields accession,organism-name,organism-tax-id,assmstats-total-sequence-len,assmstats-total-number-of-chromosomes \
			  > genomes/db_accession.tsv
	    """
}

process GENOMES_TO_FASTA {
	input:
		path('genomes')
	script:
	"""
	./db_build.R && chmod -R a+r db
	"""
}


workflow ORGFINDER_DB_DOWNLOAD {
	main:
		def taxdump = NCBI_TAXDUMP_DOWNLOAD()
		def genomes_ch = Channel.of(
			"taxon 'Pseudomonas aeruginosa'  --reference --assembly-level complete",
			"taxon 'Acinetobacter baumannii' --reference --assembly-level complete",
			"taxon 'Enterococcus'            --reference --assembly-level complete",
			"taxon 'Staphylococcus'          --reference --assembly-level complete",
			"taxon 'Streptococcus'           --reference --assembly-level complete",
			"taxon 'Enterobacterales'        --reference --assembly-level complete",
			"taxon 'Aeromonas'               --reference --assembly-level complete",
			"taxon 'Enterococcus faecalis'   --reference",
			"taxon 'Citrobacter murliniae'   --reference"
		)
		| NCBI_DATASET_DOWNLOAD_GENOME
		genomes_ch = GENOMES_AGGREGATE(genomes_ch.collect()).map({["all_collected_genomes",it]})
		RSCRIPT(genomes_ch,file("${moduleDir}/assets/db_build.R"),taxdump.view())

	emit:
		fa = Channel.empty()
}

