

include { FASTANI } from './modules/fastani'


process NCBI_DATASET_DOWNLOAD_GENOME {
    container 'docker.io/staphb/ncbi-datasets:18.18.0'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
    	val args
    output:
  		tuple path('dataset')
    script:
	    """
	    datasets download genome ${args} --filename dataset.zip \
	    	&& unzip dataset.zip -d dataset
	    """
}


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
		
		GENOMES_AGGREGATE(genomes_ch.collect())
		
	emit:
		fa = Channel.empty()
}


/*
workflow ORGFINDER {
	take:
		fa_ch
	main:
		ORGFINDER_DB_DOWNLOAD()
		FASTANI(fa_ch)
		
	emit:

}
*/

