
include { NCBI_DATASET_DOWNLOAD_GENOME } from './modules/ncbi/dataset/main.nf'
include { NCBI_TAXDUMP_DOWNLOAD        } from './modules/ncbi/taxdump/main.nf'
include { RSCRIPT                      } from './modules/rscript/main.nf'
include { FASTANI } from './modules/fastani/main.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


// Make the tsv file with all accession numbers
process ACCESSIONS_TSV {
    container 'docker.io/staphb/ncbi-datasets:18.18.0'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
  		path('genomes/query*')
  	output:
  		tuple(path('genomes/',type: 'dir'), path('accessions.tsv'))
    script:
	    """
			cat genomes/query*/ncbi_dataset/data/assembly_data_report.jsonl \
			  | dataformat tsv genome --force --fields accession,organism-name,organism-tax-id,assmstats-total-sequence-len,assmstats-total-number-of-chromosomes \
			  | sort \
			  | uniq \
			  > accessions.tsv
	    """
}

process MOLECULE_TYPES_TSV {
    container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
    memory '2 GB'
    cpus 1
    time '30 min'
    input:
  		tuple(path('ncbi_db/genomes'),path('ncbi_db/accessions.tsv'))
  	output:
  		path('ncbi_db/',type: 'dir')
    script:
	    """
			jq -r '[.assemblyAccession, .chrName, .assignedMoleculeLocationType] | @tsv' ncbi_db/genomes/query*/ncbi_dataset/data/*/sequence_report.jsonl \
			| sort \
			| uniq \
			> ncbi_db/molecule_types.tsv
	    """
}

workflow ORGFINDER_DB {
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
			"taxon 'Citrobacter murliniae'   --reference --include genome,seq-report",
			'accession GCA_040096145.1', //E. intestinihominis	3133180	GCA_040096145.1	Genome of type strain CLA-AC-H004ᵀ
			'accession GCA_001875655.1', //E. hormaechei subsp. hormaechei	301105 / parent 158836	GCA_001875655.1	Genome of species type strain ATCC 49162ᵀ
			'accession GCA_001729745.1', //E. hormaechei subsp. hoffmannii	1812934	GCA_001729745.1	Complete genome of type strain DSM 14563ᵀ
			'accession GCF_048568405.1', //E. intestinihominis RefSeq representative	3133180	GCF_048568405.1	Our current reference
			'accession GCF_019048245.1'  //E. hormaechei RefSeq representative	158836	GCF_019048245.1	Our current reference
		)
		| NCBI_DATASET_DOWNLOAD_GENOME
		
		def db_ch = genomes_ch
			.collect()
			| ACCESSIONS_TSV
			| MOLECULE_TYPES_TSV

		RSCRIPT(db_ch.map({["all",it]}),file("${moduleDir}/assets/db_build.R"),taxdump)
	emit:
		db = RSCRIPT.out.map({it[1]})
}


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




