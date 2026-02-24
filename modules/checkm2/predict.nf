
process CHECKM2_PREDICT {
	  container 'community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d'
    memory '10 GB'
    cpus 4
    time '30 min'
    input:
		    tuple val(meta), path('assembly.fasta')
		    path('db')
    output:
		    tuple val(meta), path("output/quality_report.tsv")
		script:
	    """
	    checkm2 predict \\
	        --input 'assembly.fasta' \\
	        --output-directory 'output' \\
	        --threads ${task.cpus} \\
	        --database_path ${db}/CheckM2_database/uniref100.KO.1.dmnd \\
	        ${task.ext.args?:''}
    """

    stub:
    """
    mkdir -p output/ && touch output/quality_report.tsv
    """
}

