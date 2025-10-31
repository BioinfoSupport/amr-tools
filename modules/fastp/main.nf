process FASTP {
    container 'community.wave.seqera.io/library/fastp:1.0.1--c8b87fe62dcc103c'
    memory '4 GB'
    cpus 4
    time '30 min'
    input:
	    tuple val(meta), path(reads)
    output:
	    tuple val(meta), path('*.json')           , emit: json
	    tuple val(meta), path('*.html')           , emit: html
    script:
    	if (reads.size()==1) {
		    """
		    fastp \\
		        --in1 ${reads[0]} \\
		        --thread ${task.cpus} \\
		        --json output.fastp.json \\
		        --html output.fastp.html \\
		        ${task.ext.args?:''}
		    """
    	} else {
		    """
		    fastp \\
		        --in1 ${reads[0]} --in2 ${reads[1]} \\
		        --thread ${task.cpus} \\
		        --json output.fastp.json \\
		        --html output.fastp.html \\
		        ${task.ext.args?:''}
		    """    		
    	}
}