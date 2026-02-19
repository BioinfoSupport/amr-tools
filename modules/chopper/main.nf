process CHOPPER {
    container 'quay.io/biocontainers/chopper:0.12.0--hcdda2d0_0'
    memory '8 GB'
    cpus 4
    time '30 min'
    input:
	    tuple val(meta), path('reads.fastq.gz')
    output:
	    tuple val(meta), path('chopped_reads.fastq.gz')
    script:
	    """
	    chopper ${task.ext.args?:'--quality 15 --minlength 200'} \\
	        --input reads.fastq.gz \\
	        --threads ${task.cpus} \\
	    | gzip --fast \\
	    > chopped_reads.fastq.gz
	    """
    stub:
	    """
	    echo "" | gzip > chopped_reads.fastq.gz
	    """
}
