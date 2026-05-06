

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


workflow ORGFINDER_DB_DOWNLOAD {
	main:
		def taxdump = NCBI_TAXDUMP_DOWNLOAD()
		def genomes_ch = Channel.of(
			"taxon 'Pseudomonas aeruginosa'  --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Acinetobacter baumannii' --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Enterococcus'            --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Staphylococcus'          --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Streptococcus'           --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Enterobacterales'        --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Aeromonas'               --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Myroides'                --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Enterococcus faecalis'   --reference --include genome,seq-report",
			"taxon 'Citrobacter murliniae'   --reference --include genome,seq-report"
		)
		| NCBI_DATASET_DOWNLOAD_GENOME
		
		// Filter only chromosomic sequences
		//jq -r 'select(.assignedMoleculeLocationType == "Chromosome") | (.assemblyAccession + ":" + .chrName)' work/55/e3f8291f4f02516b545ef7239acce4/dataset/ncbi_dataset/data/*/sequence_report.jsonl

		
		genomes_ch = GENOMES_AGGREGATE(genomes_ch.collect()).map({["all_collected_genomes",it]})
		RSCRIPT(genomes_ch,file("${moduleDir}/assets/db_build.R"),taxdump)
	emit:
		db = RSCRIPT.out.map({it[1]})
}

