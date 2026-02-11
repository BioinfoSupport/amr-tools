
process SAMTOOLS_FAIDX {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 1
    time '15 min'
    input:
    		tuple val(meta), path(fasta)
    output:
    		tuple val(meta), path ("*.fai")
    script:
		    """
		    samtools faidx ${task.ext.args?:''} '${fasta}'
		    """
		stub:
				"""
				touch '${fasta}.fai'
				"""
}
