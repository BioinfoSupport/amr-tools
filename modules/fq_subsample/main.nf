
process FQ_SUBSAMPLE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 4
    time '30 min'
    input:
    		tuple val(meta), path(reads)
    output:
    		tuple val(meta), path("subsampled_reads_*.fastq.gz")
    script:
    		reads = reads instanceof List?reads:[reads]
    		def bp_limit = task.ext.bp_limit?:-1
    		if (bp_limit<0) {
    			reads.withIndex().collect({x,i -> "ln -s '${x}' 'subsampled_reads_${i+1}.fastq.gz'"}).join("\n")
    		} else {
    			//TODO: When pair-end reads are given, we should sum reads length of all pairs and stop when limit is reached
    			reads.withIndex().collect({x,i -> """
							# Set this option because awk can exit and break the pipe
							set +o pipefail
							bgzip -dc '${x}' \\
							| awk 'NR%4==2{bp+=length(\$0)} {print} ((bp>=${bp_limit}) && (NR%4==0)){exit 0}' \\
							| bgzip -@ ${task.cpus} \\
							> 'subsampled_reads_${i+1}.fastq.gz'
    			"""}).join("\n")
    		}
    stub:
    	reads.withIndex().collect({x,i -> "touch 'subsampled_reads_${i+1}.fastq.gz'"}).join("\n")
}

