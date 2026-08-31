

include { NCBI_DATASET_DOWNLOAD_GENOME } from './modules/ncbi/dataset/main.nf'
include { NCBI_TAXDUMP_DOWNLOAD        } from './modules/ncbi/taxdump/main.nf'
include { RSCRIPT                      } from './modules/rscript/main.nf'


// Make the tsv file with all accession numbers
process MAKE_DB_ACCESSION_TSV {
    container 'docker.io/staphb/ncbi-datasets:18.18.0'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
  		path('genomes/query*')
  	output:
  		path('genomes/')
    script:
	    """
			cat genomes/query*/ncbi_dataset/data/assembly_data_report.jsonl \
			  | dataformat tsv genome --force --fields accession,organism-name,organism-tax-id,assmstats-total-sequence-len,assmstats-total-number-of-chromosomes \
			  > genomes/db_accession.tsv
	    """
}


process MAKE_DB_MOLECULE_TYPE_TSV {
    container 'community.wave.seqera.io/library/fastp:1.0.1--c8b87fe62dcc103c'
    memory '2 GB'
    cpus 1
    time '30 min'
    input:
  		path('genomes/query*')
  	output:
  		path('genomes/')
    script:
	    """
			jq -r '(.assemblyAccession + ":" + .chrName + "" + .assignedMoleculeLocationType)' genomes/query*/ncbi_dataset/data/*/sequence_report.jsonl \
			> genomes/assigned_molecule_types.tsv
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
		
		genomes_ch = genomes_ch.collect()
		| MAKE_DB_ACCESSION_TSV
		| MAKE_DB_MOLECULE_TYPE_TSV
		| map({["all_collected_genomes",it]})

		
		RSCRIPT(genomes_ch,file("${moduleDir}/assets/db_build.R"),taxdump)
	emit:
		db = RSCRIPT.out.map({it[1]})
}

