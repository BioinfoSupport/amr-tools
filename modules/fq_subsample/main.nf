
process FQ_SUBSAMPLE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '6 GB'
    cpus 4
    time '30 min'
    input:
    		tuple val(meta), path("reads.fastq.gz")
    		val(bp_limit)
    output:
    		tuple val(meta), path("subsampled_reads.fastq.gz")
    script:
				"""
				gzip -dc reads.fastq.gz \\
				| awk 'NR%4==2{bp+=length(\$0)} {print} ((bp>=${bp_limit}) && (NR%4==0)){exit}' \\
				| bgzip -@ ${task.cpus} \\
				> subsampled_reads.fastq.gz
				"""
    stub:
	    """
	    touch subsampled_reads.fastq.gz
	    """
}


