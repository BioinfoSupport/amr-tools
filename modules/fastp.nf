process FASTP {
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/527b18847a97451091dba07a886b24f17f742a861f9f6c9a6bfb79d4f1f3bf9d/data'
    memory '4 GB'
    cpus 4
    time '30 min'
    input:
	    tuple val(meta), path(reads)
    output:
	    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
	    tuple val(meta), path('*.json')           , emit: json
	    tuple val(meta), path('*.html')           , emit: html
	    tuple val(meta), path('*.log')            , emit: log
	    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
	    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
	    path "versions.yml"                       , emit: versions
    script:
	    """
	    fastp \\
	        --in1 read_1.fastq.gz \\
	        --in2 read_2.fastq.gz \\
	        --out1 output_1.fastp.fastq.gz \
	        --out2 output_2.fastp.fastq.gz \
	        --thread ${task.cpus} \\
	        --json output.fastp.json \\
	        --html output.fastp.html \\
	        --detect_adapter_for_pe \\
	        --failed_out output.paired.fail.fastq.gz
	        ${task.ext.args?:''} \\
	        2>| >(tee output.fastp.log >&2)
	    """
}