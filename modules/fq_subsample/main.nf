
process FQ_SUBSAMPLE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 4
    time '15 min'
    ext.limit = 1000000000
    input:
    		tuple val(meta), path("reads.fastq.gz")
    output:
    		tuple val(meta), path("subsampled_reads.fastq.gz")
    script:
				"""
				gzip -dc reads.fastq.gz \\
				| awk 'NR%4==2{bp+=length(\$0)} {print} ((bp>=${task.ext.limit}) && (NR%4==0)){exit}' \\
				| bgzip -@ ${task.cpus} \\
				> subsampled_reads.fastq.gz
				"""
}


